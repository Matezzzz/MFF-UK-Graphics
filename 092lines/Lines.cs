using System;
using System.Collections.Generic;
using System.Drawing;
using LineCanvas;
using Utilities;
using MathSupport;


namespace _092lines
{

  class GloriousDraw
  {
    static double percent (int step, int step_count) => 1.0 * step / step_count;
    static double angle (int step, int step_count, int shift=0) => 2 * Math.PI * percent(step + shift, step_count);
    static double color (int step, int step_count, int shift) => 360 * percent(step + shift, step_count);


    static void drawLine (Canvas c, int x, int y, double a1, double a2, double r1, double r2)
    {
      c.Line(x + Math.Sin(a1) * r1, y + Math.Cos(a1) * r1, x + Math.Sin(a2) * r2, y + Math.Cos(a2) * r2);
    }
    public static void DrawCircle (Canvas c, int x, int y, int steps, double r1, double r2, int color_hue_shift, int angle_shift_steps, double color_lightness)
    {
      for (int k = 0; k < steps; k++)
      {
        c.SetColor(Arith.HSVToColor(color(k, steps, color_hue_shift), 1.0, color_lightness));
        double a = angle(k, steps);
        drawLine(c, x, y, a, angle(k, steps, angle_shift_steps), r1, r2);
        drawLine(c, x, y, a, angle(k, steps,-angle_shift_steps), r1, r2);
      }
    }

    public static void DrawSpiral (Canvas c, int x, int y, int steps, double rdif, int color_hue_shift, int angle_shift_steps, int init_angle_steps, double color_lightness)
    {
      for (int k = 0; k < 10 * steps; k++)
      {
        c.SetColor(Arith.HSVToColor(color(k, steps, color_hue_shift), 1.0, color_lightness));
        int ang_k = init_angle_steps + k;
        double a = angle(ang_k, steps);
        drawLine(c, x, y, a, angle(ang_k, steps, angle_shift_steps), k + rdif, k);
        drawLine(c, x, y, a, angle(ang_k, steps,-angle_shift_steps), k + rdif, k);
      }
    }
  }

  public class Lines
  {
    /// <summary>
    /// Form data initialization.
    /// </summary>
    /// <param name="name">Your first-name and last-name.</param>
    /// <param name="wid">Initial image width in pixels.</param>
    /// <param name="hei">Initial image height in pixels.</param>
    /// <param name="param">Optional text to initialize the form's text-field.</param>
    /// <param name="tooltip">Optional tooltip = param help.</param>
    public static void InitParams (out string name, out int wid, out int hei, out string param, out string tooltip)
    {
      // {{

      // Put your name here.
      name = "Matěj Mrázek";

      // Image size in pixels.
      wid = 2560;
      hei = 1440;

      // Specific animation params.
      param = "width=1.0,anti=true";

      // Tooltip = help.
      tooltip = "width=<int>, anti[=<bool>]";

      // }}
    }




    



    /// <summary>
    /// Draw the image into the initialized Canvas object.
    /// </summary>
    /// <param name="c">Canvas ready for your drawing.</param>
    /// <param name="param">Optional string parameter from the form.</param>
    public static void Draw (Canvas c, string param)
    {
      // {{ TODO: put your drawing code here

      // Input params.
      float penWidth = 1.0f;   // pen width
      bool antialias = false;  // use anti-aliasing?

      Dictionary<string, string> p = Util.ParseKeyValueList(param);
      if (p.Count > 0)
      {
        // with=<line-width>
        if (Util.TryParse(p, "width", ref penWidth))
        {
          if (penWidth < 0.0f)
            penWidth = 0.0f;
        }

        // anti[=<bool>]
        Util.TryParse(p, "anti", ref antialias);

      }

      int c_wid = c.Width;
      int c_hei = c.Height;
      int c_siz = Math.Min(c_wid, c_hei);

      c.Clear(Color.Black);
      c.SetAntiAlias(antialias);
      c.SetPenWidth(penWidth);


      int largest_circle_r = (c_siz - 30);

      for (float i = largest_circle_r; i > 10; i /= 1.1F)
      {
        GloriousDraw.DrawCircle(c, c_wid / 2, c_hei / 2, 72, i, i / 1.1F, (int) i, 3, 0.25F);
      }
      for (float i = 0; i < 4; i++)
      {
        GloriousDraw.DrawSpiral(c, c_wid / 2, c_hei / 2, 72, 10, (int) i + 100, 3, 18 * (int) i, 0.4F);
      }
      for (float i = 20; i < largest_circle_r * 0.6; i *= 1.8F)
      {
        GloriousDraw.DrawCircle(c, c_wid / 2, c_hei / 2, 360, i, i / 1.2, 0, 45, 0.6);
      }
    }
  }
}
