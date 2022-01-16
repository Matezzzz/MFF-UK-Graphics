using MathSupport;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Globalization;
using System.Threading.Tasks;
using Utilities;

namespace Modules
{
  public class ModulePaint : DefaultRasterModule
  {
    /// <summary>
    /// Mandatory plain constructor.
    /// </summary>
    public ModulePaint ()
    {
      // Default HSV transform parameters (have to by in sync with the default module state).
      param = "";
    }

    /// <summary>
    /// Author's full name.
    /// </summary>
    public override string Author => "Matěj Mrázek";

    /// <summary>
    /// Name of the module (short enough to fit inside a list-boxes, etc.).
    /// </summary>
    public override string Name => "Paint filter";

    /// <summary>
    /// Tooltip for Param (text parameters).
    /// </summary>
    public override string Tooltip =>
      "[,gamma=<double>][,slow][,par]";

    /// <summary>
    /// Usually read-only, optionally writable (client is defining number of inputs).
    /// </summary>
    public override int InputSlots => 1;

    /// <summary>
    /// Usually read-only, optionally writable (client is defining number of outputs).
    /// </summary>
    public override int OutputSlots => 1;

    /// <summary>
    /// Input raster image.
    /// </summary>
    protected CoolImage in_image;

    /// <summary>
    /// Output raster image.
    /// </summary>
    protected CoolImage out_image;
    protected CoolImage grad_image;
    protected CoolImage grad2_image;
    protected CoolImage grad_blurred;

    /// <summary>
    /// Slow computation (using GetPixel/SetPixel).
    /// </summary>
    public bool slow = false;

    /// <summary>
    /// Parallel computation (most useful for large pictures).
    /// </summary>
    public bool parallel = true;

    protected void updateParam ()
    {
      if (paramDirty)
      {
        paramDirty = false;

        // 'param' parsing.
        Dictionary<string, string> p = Util.ParseKeyValueList(param);
        if (p.Count > 0)
        {
          // dH=<double> [+- number in degrees]
          //Util.TryParse(p, "dH", ref dH);

         

          // par .. use Parallel.For
          parallel = p.ContainsKey("par");

          // slow .. set GetPixel/SetPixel computation
          slow = p.ContainsKey("slow");
        }
      }
    }


    protected class CoolImage
    {
      public Bitmap bitmap;
      BitmapData data;
      int stride, pixel_bytes;
      public CoolImage (int w, int h) : this(new Bitmap(w, h, PixelFormat.Format24bppRgb))
      {
      }
      public CoolImage(Bitmap b)
      {
        bitmap = b;
      }
      public int Width
      {
        get => bitmap.Width;
      }
      public int Height
      {
        get => bitmap.Height;
      }
      public void Lock (ImageLockMode mode)
      {
        data = bitmap.LockBits(new Rectangle(0, 0, bitmap.Width, bitmap.Height), mode, bitmap.PixelFormat);
        stride = data.Stride;
        pixel_bytes = Image.GetPixelFormatSize(data.PixelFormat) / 8;
      }
      public void LockR ()
      {
        Lock(ImageLockMode.ReadOnly);
      }
      public void LockW ()
      {
        Lock(ImageLockMode.WriteOnly);
      }
      public void LockRW ()
      {
        Lock(ImageLockMode.ReadWrite);
      }
      public void Unlock ()
      {
        bitmap.UnlockBits(data);
      }
      unsafe public byte* getPtr (int x, int y)
      {
        return (byte*)(data.Scan0 + x * pixel_bytes + y * stride);
      }
    }


    class DColor
    {
      //RGB are doubles on the interval <0, 1)
      public double[] RGB;
      public double R
      {
        get => RGB[0];
        set => RGB[0] = value;
      }
      public double G
      {
        get => RGB[1];
        set => RGB[1] = value;
      }
      public double B
      {
        get => RGB[2];
        set => RGB[2] = value;
      }
      public double Saturation => Math.Max(R, Math.Max(G, B)) - Math.Min(R, Math.Min(G, B));

      public DColor (DColor c) => RGB = new double[3] { c.R, c.G, c.B };

      public DColor (double r, double g, double b) => RGB = new double[3] { r, g, b };

      public void MakeAbs ()
      {
        R = Math.Abs(R);
        G = Math.Abs(G);
        B = Math.Abs(B);
      }

      unsafe public DColor (byte* vals)
      {
        RGB = new double[3];
        Set(vals);
      }

      unsafe public void Set (byte* vals)
      {
        R = vals[2] / 255.0;
        G = vals[1] / 255.0;
        B = vals[0] / 255.0;
      }

      unsafe public void SaveTo (byte* ptr)
      {
        ptr[2] = Convert.ToByte(Clamp(R) * 255);
        ptr[1] = Convert.ToByte(Clamp(G) * 255);
        ptr[0] = Convert.ToByte(Clamp(B) * 255);
      }


      public static DColor FromColor (Color c)
      {
        return new DColor(c.R / 256.0, c.G / 256.0, c.B / 256.0);
      }
      public static DColor operator + (DColor c1, DColor c2) => new DColor(c1.R + c2.R, c1.G + c2.G, c1.B + c2.B);
      public static DColor operator - (DColor c1, DColor c2) => new DColor(c1.R - c2.R, c1.G - c2.G, c1.B - c2.B);
      public static DColor operator * (DColor c1, double x) => new DColor(c1.R * x, c1.G * x, c1.B * x);
      public static DColor operator / (DColor c1, double x) => new DColor(c1.R / x, c1.G / x, c1.B / x);
      public Color ToColor ()
      {
        return Color.FromArgb((int)(Clamp(R) * 255), (int)(Clamp(G) * 255), (int)(Clamp(B) * 255));
      }

      public void reset ()
      {
        R = 0;
        G = 0;
        B = 0;
      }



      //think of the color as a vector and compute its' magnitude
      public double VecMagnitude ()
      {
        return System.Math.Sqrt(VecMagnitudeSquared());
      }

      //think of the color as a vector and compute its' magnitude squared
      public double VecMagnitudeSquared ()
      {
        return R * R + G * G + B * B;
      }

      //clamp value between 0 and 1
      private static double Clamp (double color)
      {
        return color < 0 ? 0 : (color > 1 ? 1 : color);
      }

      static double r_, dR, dG, dB, diff;
      //distance between two colors. Should be between 0 and 1. Equation taken from wikipedia https://en.wikipedia.org/wiki/Color_difference
      public double RedmeanDistanceFrom (DColor c)
      {
        r_ = (R + c.R) / 2;
        dR = R - c.R;
        dG = G - c.G;
        dB = B - c.B;
        diff = ((2 + r_) * dR * dR + 4 * dG * dG + (3.0 - r_) * dB * dB) / 9.0;
        return Math.Sqrt(diff);
      }

      //compute distance between two colors, euclidean
      public double EuclideanDistanceFrom (DColor c)
      {
        dR = R - c.R;
        dG = G - c.G;
        dB = B - c.B;
        diff = (dR * dR + dG * dG + dB * dB) / 3;
        return Math.Sqrt(diff);
      }

      //compute distance between two colors based on sum of absolute values of differences
      public double AbsoluteDistanceFrom (DColor c)
      {
        dR = R - c.R;
        dG = G - c.G;
        dB = B - c.B;
        return (Math.Abs(dR) + Math.Abs(dG) + Math.Abs(dB)) / 3;
      }

      //this checks for equality of all components
      //as these are doubles, this method is only usable for colors that were created the same way(like from the same uint8 RGB color)
      //in most other cases, this will simply be false due to numerical errors
      public static bool operator == (DColor c1, DColor c2) => c1.R == c2.R && c1.G == c2.G && c1.B == c2.B;
      //inverse of the one above
      public static bool operator != (DColor c1, DColor c2) => !(c1 == c2);


      //if this method isn't overriden, compiler complains, so here it is
      public override bool Equals (object obj)
      {
        return base.Equals(obj);
      }

      //if this method isn't overriden, compiler complains, so here it is
      public override int GetHashCode ()
      {
        return base.GetHashCode();
      }
    }
  

    /// <summary>
    /// Assigns an input raster image to the given slot.
    /// Doesn't start computation (see #Update for this).
    /// </summary>
    /// <param name="inputImage">Input raster image (can be null).</param>
    /// <param name="slot">Slot number from 0 to InputSlots-1.</param>
    public override void SetInput (
      Bitmap inputImage,
      int slot = 0)
    {
      in_image = new CoolImage(inputImage);
    }

    /// <summary>
    /// Recompute the output image[s] according to input image[s].
    /// Blocking (synchronous) function.
    /// #GetOutput() functions can be called after that.
    /// </summary>
    public override void Update ()
    {
      if (in_image == null)
        return;

      // Update module values from 'param' string.
      updateParam();

      int wid = in_image.Width;
      int hei = in_image.Height;

      // Output image must be true-color.
      out_image = new CoolImage(wid, hei);
      grad_image = new CoolImage(wid, hei);
      grad2_image = new CoolImage(wid, hei);
      grad_blurred = new CoolImage(wid, hei);

      if (slow)
      {
        // Slow GetPixel/SetPixel code.

        throw new NotImplementedException("Sorry, we don't do slow.");
      }
      else
      {
        // Fast memory-mapped code.
        PixelFormat iFormat = in_image.bitmap.PixelFormat;
        if (!PixelFormat.Format24bppRgb.Equals(iFormat) &&
            !PixelFormat.Format32bppArgb.Equals(iFormat) &&
            !PixelFormat.Format32bppPArgb.Equals(iFormat) &&
            !PixelFormat.Format32bppRgb.Equals(iFormat))
          iFormat = PixelFormat.Format24bppRgb;

        Rectangle rect = new Rectangle(0, 0, wid, hei);
        in_image.LockR();
        out_image.LockW();
        grad_image.LockRW();
        grad2_image.LockRW();
        grad_blurred.LockRW();
        
        unsafe
        {
          int dI = Image.GetPixelFormatSize(iFormat) / 8;

          

          void KernelMethod(int y, double[,] kernel, CoolImage data_in, CoolImage data_out, bool apply_abs = false)
          {
            DColor cur_color = new DColor(0, 0, 0);
            DColor mid_color = new DColor(0, 0, 0);
            DColor shift = new DColor(0.5, 0.5, 0.5);
            DColor result = new DColor(0, 0, 0);
            int wx, wy;

            int from_x = -kernel.GetLength(0) / 2;
            int from_y = -kernel.GetLength(1) / 2;

            for (int x = 0; x < wid; x++)
            {
              byte* ptr = data_in.getPtr(x, y);
              mid_color.Set(ptr);
              result.reset();
              for (int kx = 0; kx < kernel.GetLength(0); kx++)
              {
                for (int ky = 0; ky < kernel.GetLength(1); ky++)
                {
                  wx = x + kx + from_x;
                  wy = y + ky + from_y;
                  if (0 > wx || wx >= wid || 0 > wy || wy >= hei)
                    continue;
                  ptr = data_in.getPtr(wx, wy);
                  cur_color.Set(ptr);
                  result += cur_color * kernel[kx,ky];
                }
              }
              if (apply_abs)
              {
                result.MakeAbs();
              }

              ptr = data_out.getPtr(x, y);
              result.SaveTo(ptr);
            }
          }


          double[,] grad_kernel = new double[3,3]
          {
            {1, 2, -1},
            {2, 0, -2},
            {1, -2, -1}
          };
          void gradient(int y, CoolImage data_in, CoolImage grad_out)
          {
            KernelMethod(y, grad_kernel, data_in, grad_out, true);
          }



          void gradient1(int y)
          {
            gradient(y, in_image, grad_image);
          }
          void gradient2(int y)
          {
            gradient(y, grad_image, grad2_image);
          }


          double[,] blur_kernel = new double[3,3]
          {
            {0.0625, 0.125, 0.0625},
            {0.125, 0.25, 0.125},
            {0.0625, 0.125, 0.0625}
          };
          void blur(int y, CoolImage data_in, CoolImage data_out)
          {
            KernelMethod(y, blur_kernel, data_in, data_out);
          }

          void blurGradient(int y)
          {
            blur(y, grad2_image, grad_blurred);
          }


          void inner (int y)
          {
            // User break handling.
            if (UserBreak)
              return;

            double S1, S2;
            double H, V1, V2;
            byte* ptr;
            DColor cur = new DColor(0, 0, 0);
            DColor result;
            DColor dcolor_black = new DColor(0, 0, 0);
            for (int x = 0; x < wid; x++)
            {
              // Recompute one pixel (*iptr -> *optr).
              // iptr, optr -> [B,G,R]

              DColor current = new DColor(in_image.getPtr(x, y));
              DColor grad = new DColor(grad_blurred.getPtr(x, y));
              DColor grad2 = new DColor(grad2_image.getPtr(x, y));

              int compareColors (DColor a, DColor b)
              {
                Arith.RGBtoHSV(a.R, a.G, a.B, out H, out S1, out V1);
                Arith.RGBtoHSV(b.R, b.G, b.B, out H, out S2, out V2);
                return V1.CompareTo(V2);
              }


              DColor Median (List<DColor> colors)
              {
                colors.Sort(compareColors);
                return colors[(colors.Count - 1) / 2];
              }

              List <DColor> l = new List<DColor>();

              int width_half = 3;
              int height_half = 3;

              int wx, wy;
              for (int dx = -width_half; dx <= width_half; dx++)
              {
                for (int dy = -height_half; dy <= height_half; dy++)
                {
                  wx = x + dx; wy = y + dy;
                  if (0 > wx || wx >= wid || 0 > wy || wy >= hei) continue;
                  ptr = in_image.getPtr(wx, wy);
                  DColor color = new DColor(ptr);
                  l.Add(new DColor(ptr));
                }
              }

              double grad2_mag = grad2.VecMagnitude() / Math.Sqrt(3);
              double grad_mag = grad.VecMagnitude() / Math.Sqrt(3) * grad2_mag;
              double grad_coef = grad_mag;
              if (grad_mag > 0.25)
              {
                double c = Math.Min(1, (grad_coef - 0.25) / 0.25);
                result = Median(l) * (1.0 - c) + dcolor_black * c;
              }
              else
              {
                result = Median(l);
              }
              result.SaveTo(out_image.getPtr(x, y));
            }
            
          }
          void copy (int y, CoolImage src, CoolImage dst)
          {
            for (int x = 0; x < wid; x++)
            {
              byte* s = src.getPtr(x, y);
              byte* d = dst.getPtr(x, y);
              d[0] = s[0];
              d[1] = s[1];
              d[2] = s[2];
            }
          }

          void save(int y, CoolImage src)
          {
            copy(y, src, out_image);
          }

         void saveGradient2 (int y)
          {
            save(y, grad_blurred);
          }



          if (parallel)
          {
            _=Parallel.For(0, hei, gradient1);
            _=Parallel.For(0, hei, gradient2);
            _ = Parallel.For(0, hei, blurGradient);

            //_=Parallel.For(0, hei, inner);
            _=Parallel.For(0, hei, saveGradient2);
          }
          else
          {
            for (int y = 0; y < hei; y++)
            {
              inner(y);
            }
          }
        }

        in_image.Unlock();
        out_image.Unlock();
        grad_image.Unlock();
        grad2_image.Unlock();
        grad_blurred.Unlock();
      }
    }
    class InnerFunc
    {

    }

    /// <summary>
    /// Returns an output raster image.
    /// Can return null.
    /// </summary>
    /// <param name="slot">Slot number from 0 to OutputSlots-1.</param>
    public override Bitmap GetOutput (int slot = 0) => out_image.bitmap;
  }
}
