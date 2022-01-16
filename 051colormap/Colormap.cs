using System.Drawing;
using MathSupport;
using System.Linq;
using System;

/**
 * Main idea:
 *  1.  Try many random combinations of colors, start with best ones (the ones with the lowest error function)
 *  2.  Loop over:
 *  3.    Compute matrix of gradients of error function. The input to this function is the current palette of colors.
 *          The matrix looks as follows - let d denote the partial derivative, E the error function, Ri...Bi the components of the i-th palette color:
 *            (dE / dR1, dE / dG1, dE / dB1)
 *            (dE / dR2, dE / dG2, dE / dB2)
 *            (   .         .         .    )
 *            (dE / DRn, dE / dGn, dE / dBn)
 *          n is number of palette colors.
 *          The error function can take on many forms (distance from a palette means distance to the closest palette color):
 *            It can be exact - sum of distances of all image pixels, divided by the pixel count (=avg. distance for sample)
 *            It can be sampled - sum of distances from N samples distributed at random over the image (=avg.distance for sample)
 *            It can be computed from a color cube - create a cube, K by K by K colors, each point holds count of pixels closest to it.
 *              During computing distance, go through all the points, compute distance each one, then compute the weighted average using pixel counts for each point.
 *          Distance between colors used in the error function above can either be computed using the euclidean distance function or using the redmean distance, which accounts for how humans view colors. 
 *  5.    Perform gradient descent using the matrix described above - subtract each gradient from its' respective color
 *  6.    If there is a color with which zero samples are associated, rewrite it with a new random color. This is done to eliminate useless colors in the final palette.      
 *
 *  Also, if there are less colors in the image than required in the palette, just pick those instead of running the algorithm above.
 */



namespace _051colormap
{

  class Colormap
  {
    /// <summary>
    /// Form data initialization.
    /// </summary>
    public static void InitForm (out string author)
    {
      author = "Matěj Mrázek";
    }


    //if true, convert the image to match the palette after it is computed. Wasn't in the task, but was easy to make, and it looks cool
    const bool convert_image_to_palette = false;
    //taking all pixels into account when computing distance is too slow, so N samples are taken, and distance is measured from those
    const int distance_sample_count = 1000;
    
    enum DistanceFunction
    {
      EXACT, SAMPLES, COLOR_CUBE
    };
    //which distance function will be used to compute distance from the image
    const DistanceFunction distance_function = DistanceFunction.COLOR_CUBE;

    //how many subidivisions are used for each color during the color cube distance measuring method
    const int color_cube_subdivisions = 16;


    enum ColorDistanceFunction
    {
      REDMEAN, EUCLIDEAN, ABSOLUTE
    };
    //Redmean takes into account how humans view colors. Euclidean is normalized distance between two vectors. Absolute is the normalized sum of absolute values of differences.
    const ColorDistanceFunction color_distance_function = ColorDistanceFunction.ABSOLUTE;

    //how many palettes are solved for before the best one is selected
    const int palette_test_count = 1;



    //class that measures how well a given palette matches an image
    class ImageDistance
    {
      //image to compute distance from
      readonly DImage image;

      //positions of samples
      readonly (int, int)[] sample_positions;

      //How many pixels are associated with each point in the color cube distance method
      readonly int[,,] color_counts;

      public ImageDistance (DImage img)
      {
        image = img;
        sample_positions = new (int, int)[distance_sample_count];
        InitSampled();
        color_counts = new int[color_cube_subdivisions, color_cube_subdivisions, color_cube_subdivisions];
        InitColorCube();
      }
      //generate new sample positions
      private void InitSampled()
      {
        for (int i = 0; i < distance_sample_count; i++)
        {
          //pick position of a random pixel in the image
          sample_positions[i] = (rand.Next(image.Width), rand.Next(image.Height));
        }
      }

      //compute how many pixels are associated with each pixel in the color cube
      private void InitColorCube ()
      {
        for (int x = 0; x < image.Width; x++)
        {
          for (int y = 0; y < image.Height; y++)
          {
            DColor c = image.Data[x, y] * color_cube_subdivisions;
            color_counts[(int)c.R, (int)c.G, (int)c.B]++;
          }
        }
      }
     

      //compute exact distance of palette from image - go through all the pixels, compute distance from each one, divide distance by pixel count
      private double DistanceFromImageExact (Palette palette)
      {
        //sum of all distances
        double total_dist = 0.0;
        for (int x = 0; x < image.Width; x++)
        {
          for (int y = 0; y < image.Height; y++)
          {
            (double, int) closest_color_dist_and_index = DistanceFromClosestPaletteColor(image.Data[x, y], palette);
            total_dist += closest_color_dist_and_index.Item1;
          }
        }
        //compute average distance per pixel
        return total_dist / image.Width / image.Height;
      }

      //compute distance of samples from the palette - go through all samples, compute distance from palette for each one, then compute their average
      public double DistanceFromImageSamples(Palette palette)
      {
        double total_dist = 0.0;
        for (int i = 0; i < distance_sample_count; i++)
        {
          int x = sample_positions[i].Item1;
          int y = sample_positions[i].Item2;
          (double, int) closest_color_dist_and_index = DistanceFromClosestPaletteColor(image.Data[x, y], palette);
          total_dist += closest_color_dist_and_index.Item1;
        }
        //compute average distance per pixel
        return total_dist / distance_sample_count;
      }

      //compute distance of palette from image using the color cube
      // = go through all cube subdivisions, select closest palette color for each one, then compute weighed average for all subdivisions based on pixel counts
      public double DistanceFromImageColorCube(Palette palette)
      {
        double total_dist = 0;
        DColor c = new DColor(0, 0, 0);
        for (int r = 0; r < color_cube_subdivisions; r++)
        {
          for (int g = 0; g < color_cube_subdivisions; g++)
          {
            for (int b = 0; b < color_cube_subdivisions; b++)
            {
              c.R = (r + 0.5) / color_cube_subdivisions;
              c.G = (g + 0.5) / color_cube_subdivisions;
              c.B = (b + 0.5) / color_cube_subdivisions;
              (double, int) closest_color_dist_and_index = DistanceFromClosestPaletteColor(c, palette);
              total_dist += closest_color_dist_and_index.Item1 * color_counts[r, g, b];
            }
          }
        }
        return total_dist / image.Width / image.Height;
      }

      //Select one of the distance functions above based on settings
      public double DistanceFromImage(Palette palette)
      {
        //mute unreachable code warning
#pragma warning disable 162
        switch (distance_function)
        {
          case DistanceFunction.EXACT:
            return DistanceFromImageExact(palette);
          case DistanceFunction.SAMPLES:
            return DistanceFromImageSamples(palette);
          case DistanceFunction.COLOR_CUBE:
            return DistanceFromImageColorCube(palette);
          default:
            Console.WriteLine("Given color function not implemented.");
        }
#pragma warning restore 162
      }

      //compute how many samples are closest to each palette color
      public int[] ImagePalleteAssociatedSampleCounts(Palette palette)
      {
        //one sample count for every palette color
        int[] sample_counts = new int[palette.ColorCount];
        for (int i = 0; i < distance_sample_count; i++)
        {
          int x = sample_positions[i].Item1;
          int y = sample_positions[i].Item2;

          (double, int) closest_color_dist_and_index = DistanceFromClosestPaletteColor(image.Data[x, y], palette);
          //increase sample count for the color this sample is closest to
          sample_counts[closest_color_dist_and_index.Item2]++;
        }
        return sample_counts;
      }

      //given color and a palette, compute the index of the closest palette color, and the distance from it
      public (double, int) DistanceFromClosestPaletteColor(DColor color, Palette palette)
      {
        double min_dist = double.PositiveInfinity;
        int min_dist_i = 0;
        for (int j = 0; j < palette.ColorCount; j++)
        {
          double dist = 0;
          //compute distance of color to current palette color - use specified distance function
#pragma warning disable 162
          switch (color_distance_function)
          {
            case ColorDistanceFunction.REDMEAN:
              dist = palette.Colors[j].RedmeanDistanceFrom(color);
              break;
            case ColorDistanceFunction.EUCLIDEAN:
              dist = palette.Colors[j].EuclideanDistanceFrom(color);
              break;
            case ColorDistanceFunction.ABSOLUTE:
              dist = palette.Colors[j].AbsoluteDistanceFrom(color);
              break;
            default:
              Console.WriteLine("Given color distance function not implemented.");
#pragma warning restore 162
          }
          if (dist < min_dist)
          {
            min_dist = dist;
            min_dist_i = j;
          }
        }
        return (min_dist, min_dist_i);
      }

      //testing method for exact distance vs. approximate distance method testing
      //ideally, distance for exact and approximate methods would be the same
      //this method computes max. error and avg. error - these decrease when increasing distance_sample_count constant
      //with default value of distance_sample_count = 1000, common max. error is around 4%, while the avg. is around 1% 
      public void DistanceApproxTest (int palette_color_count)
      {
        const int test_count = 100;

        double max_error = double.NegativeInfinity;
        double avg_error = 0.0;
        for (int i = 0; i < test_count; i++)
        {
          Palette p = Palette.GetRandom(palette_color_count);
          double real_dist = DistanceFromImageExact(p);
          double approx_dist = DistanceFromImage(p);
          double approx_accuracy = approx_dist / real_dist;
          double error = System.Math.Abs(1 - approx_accuracy);
          max_error = System.Math.Max(error, max_error);
          avg_error += error;
        }
        avg_error /= test_count;
        System.Console.WriteLine("Tests ran: {0:D}, Max error: {1:P}, Avg. error: {2:P}", test_count, max_error, avg_error);
      }

      //testing method for exact distance vs. color cube distance method testing
      //ideally, distance for exact and approximate methods would be the same
      //this method computes max. error and avg. error - these decrease when increasing distance_sample_count constant
      //with default value of color_cube_subdivisions = 8, common max. error is around 3.5%, while the avg. is around 2%
      //with ccs = 16, max. error is smaller than 1%, avg is 0.15%
      //this method is better than sampling in most cases, however, gradients can very rarely spike when being computed(point with many pixels is taken over by a new color), which makes this method a bit unstable
      public void DistanceApproxTestCube (int palette_color_count)
      {
        const int test_count = 100;

        double max_error = double.NegativeInfinity;
        double avg_error = 0.0;
        for (int i = 0; i < test_count; i++)
        {
          Palette palette = Palette.GetRandom(palette_color_count);
          double real_dist = DistanceFromImageExact(palette);
          double approx_dist = DistanceFromImageColorCube(palette);
          double approx_accuracy = approx_dist / real_dist;
          double error = Math.Abs(1 - approx_accuracy);
          max_error = Math.Max(error, max_error);
          avg_error += error;
        }
        avg_error /= test_count;
        System.Console.WriteLine("Tests ran: {0:D}, Max error: {1:P}, Avg. error: {2:P}", test_count, max_error, avg_error);
      }
    };



    /**
     * Holds all image pixels and original image reference
     * All RGB components are saved as doubles
     */
    class DImage
    {
      readonly Bitmap image;
      public int Width => Data.GetLength(0);
      public int Height => Data.GetLength(1);
  
      public DColor[,] Data{get;}

      //constructor - copy image pixel values to image_data
      public DImage (Bitmap image_)
      {
        image = image_;
        Data = new DColor[image.Width, image.Height];
        for (int x = 0; x < Width; x++)
        {
          for (int y = 0; y < Height; y++)
          {
            Data[x, y] = DColor.FromColor(image.GetPixel(x, y));
          }
        }
      }

      //return min_colors if image has more than min_colors, or the amount of colors in the image otherwise
      public int CountColorsUntil (int min_colors)
      {
        DColor[] found_colors = GetFirstNColors(min_colors);
        return found_colors.Length;
      }

      //return first N colors found in the image, searching from the top by rows
      public DColor[] GetFirstNColors(int color_count)
      {
        DColor[] found_colors = new DColor[color_count];
        int i = 0;
        for (int x = 0; x < Width; x++)
        {
          for (int y = 0; y < Height; y++)
          {
            //find out whether color is present in found colors already
            bool found = false;
            for (int j = 0; j < i; j++)
            {
              //yes, I am using == for doubles - colors have to be exactly equal to count as the same
              if (found_colors[j] == Data[x, y])
              {
                found = true;
              }
            }
            //if color isn't in found colors field yet, add it
            if (!found)
            {
              found_colors[i++] = Data[x, y];
              if (i == color_count)
              {
                return found_colors;
              }
            }
          }
        }
        //convert found colors to DColor
        DColor[] color_arr = new DColor[i];
        for (int j = 0; j < i; j++)
        {
          color_arr[j] = found_colors[j];
        }
        return color_arr;
      }

      //save modified image_data back to original bitmap
      public void SaveImageDataToBitmap ()
      {
        for (int x = 0; x < Width; x++)
        {
          for (int y = 0; y < Height; y++)
          {
            image.SetPixel(x, y, Data[x, y].ToColor());
          }
        }
      }

      //select a random pixel, then return its' color
      public DColor GetRandomPixelColor ()
      {
        return new DColor(Data[rand.Next(image.Width), rand.Next(image.Height)]);
      }
    };




    static readonly Random rand = new Random(42);

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
      public static DColor FromColor (Color c)
      {
        return new DColor(c.R / 256.0, c.G / 256.0, c.B / 256.0);
      }
      public static DColor GetRandom()
      {
        return new DColor(rand.NextDouble(), rand.NextDouble(), rand.NextDouble());
      }
      public static DColor operator + (DColor c1, DColor c2) => new DColor(c1.R + c2.R, c1.G + c2.G, c1.B + c2.B);
      public static DColor operator - (DColor c1, DColor c2) => new DColor(c1.R - c2.R, c1.G - c2.G, c1.B - c2.B);
      public static DColor operator * (DColor c1, double x) => new DColor(c1.R * x, c1.G * x, c1.B * x);
      public Color ToColor ()
      {
        return Color.FromArgb((int) (Clamp(R) * 255), (int) (Clamp(G) * 255), (int) (Clamp(B) * 255));
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
      private static double Clamp(double color)
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
      public double EuclideanDistanceFrom(DColor c)
      {
        dR = R - c.R;
        dG = G - c.G;
        dB = B - c.B;
        diff = (dR * dR + dG * dG + dB * dB) / 3;
        return Math.Sqrt(diff);
      }

      //compute distance between two colors based on sum of absolute values of differences
      public double AbsoluteDistanceFrom(DColor c)
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




    //Holds all palette colors for a given image and image distance
    class Palette
    {
      public DColor[] Colors;

      public int ColorCount => Colors.Length;

      //How many iterations of gradient descent should be used when solving for palette colors
      //I have tried plotting average sample error over time, and it rarely improves after 50 iterations, so I used this constant
      const int palette_solve_iteration_count = 50;
      //how fast does the algorithm converge, larger values can cause oscillations and instabilities
      const double gradient_descent_convergence_speed = 0.4;
      //one gradient component is equal to (E(x1, ..., xj + gradient_delta, ..., xn) - E(x1, ..., xj - gradient_delta, ..., xn)) / (2 * gradient_delta)
      //if the gradient would be computed analytically, it would be approaching zero, however, this value is used for computing a numerical approximation
      const double gradient_delta = 0.001;

      //how many different colors are tried before the best one is selected
      const int init_color_tries = 1000;

      public Palette (DColor[] colors) => Colors = colors;
      public Palette (DImage image, ImageDistance image_dist, int palette_color_count)
      {
        Colors = new DColor[palette_color_count];
        InitColors(image, image_dist);
      }
      public static Palette GetRandom (int color_count)
      {
        DColor[] colors = new DColor[color_count];
        for (int i = 0; i < color_count; i++)
        {
          colors[i] = DColor.GetRandom();
        }
        return new Palette(colors);
      }

      //put random colors from the image into the palette, try many times and select the best ones
      private void InitColors (DImage image, ImageDistance image_dist)
      {
        double best_dist = double.PositiveInfinity;
        for (int j = 0; j < init_color_tries; j++)
        {
          //select ColorCount random colors from the image
          DColor[] cur_colors = new DColor[ColorCount];
          for (int i = 0; i < ColorCount; i++)
          {
            cur_colors[i] = image.GetRandomPixelColor();
          }
          //compute distance of selected colors from image, save them if they are better than all that came before
          double d = image_dist.DistanceFromImage(new Palette(cur_colors));
          if (d < best_dist)
          {
            Colors = cur_colors;
            best_dist = d;
          }
        }
      }

      //make palette colors match the provided image
      public void SolveForPaletteColors (DImage image, ImageDistance image_dist)
      {
        for (int i = 0; i < palette_solve_iteration_count; i++)
        {
          //compute gradients for all palette colors, then subtract them from their respective colors
          PaletteSolveStep(image_dist);
          //if there is a color that has zero samples associated with it, replace it with a new random color
          FixUnusedColors(image, image_dist);
        }
        //sort colors by the amount of associated samples - most used colors at the top
        SortColors(image_dist);
      }
      private void PaletteSolveStep (ImageDistance image_dist)
      {
        //compute matrix of gradients - each element represents a gradient for each palette color
        //in each row, each element means how the error function changes with respect to each component
        DColor[] gradient = ComputeMatrixOfGradients(image_dist);
        //subtract gradient from all palette colors
        for (int j = 0; j < ColorCount; j++)
        {
          Colors[j] -= gradient[j] * gradient_descent_convergence_speed;
        }
      }

      //if palette color isn't being used, replace it with a new random color
      private void FixUnusedColors (DImage image, ImageDistance image_dist)
      {
        int[] sample_counts = image_dist.ImagePalleteAssociatedSampleCounts(this);
        for (int j = 0; j < ColorCount; j++)
        {
          //if color has no samples associated with it, rewrite it with a random color from the image, then try fixing the new set of colors
          if (sample_counts[j] == 0)
          {
            Colors[j] = image.GetRandomPixelColor();
            FixUnusedColors(image, image_dist);
            return;
          }
        }
      }

      //sort colors from the most to the least used one
      private void SortColors (ImageDistance image_dist)
      {
        int[] sample_counts = image_dist.ImagePalleteAssociatedSampleCounts(this);
        Array.Sort(sample_counts, Colors);
        Array.Reverse(Colors);
      }

      //compute gradient for each palette color and return them as a matrix
      private DColor[] ComputeMatrixOfGradients (ImageDistance image_dist)
      {
        DColor[] gradient_matrix = new DColor[ColorCount];
        for (int i = 0; i < ColorCount; i++)
        {
          gradient_matrix[i] = ComputeGradientForColor(image_dist, i);
        }
        return gradient_matrix;
      }

      //compute gradient for one palette color
      private DColor ComputeGradientForColor (ImageDistance image_dist, int color_i)
      {
        return new DColor(
          ComputeGradientForColorElement(image_dist, color_i, 0),
          ComputeGradientForColorElement(image_dist, color_i, 1),
          ComputeGradientForColorElement(image_dist, color_i, 2)
        );
      }

      //compute gradient for one palette color element
      //example equation (elem_i = 0 (R component, i-th color)): E is a function that measures a distance of palette from the image
      //  grad = (E(RGB1, ..., (Ri + gradient_delta, Gi, Bi), ..., RGBn) - E(RGB1, ..., (Ri - gradient_delta, Gi, Bi), ..., RGBn)) / (2 * gradient_delta)
      private double ComputeGradientForColorElement (ImageDistance image_dist, int color_i, int elem_i)
      {
        Colors[color_i].RGB[elem_i] += gradient_delta;
        double a = image_dist.DistanceFromImage(this);
        Colors[color_i].RGB[elem_i] -= 2 * gradient_delta;
        double b = image_dist.DistanceFromImage(this);
        Colors[color_i].RGB[elem_i] += gradient_delta;
        return (a - b) / (2 * gradient_delta);
      }

      //convert image to current palette
      public void ConvertImage(DImage image, ImageDistance image_dist)
      {
        for (int x = 0; x < image.Width; x++)
        {
          for (int y = 0; y < image.Height; y++)
          {
            //set image pixel to the closest palette color
            (double, int) closest_color_dist_and_i = image_dist.DistanceFromClosestPaletteColor(image.Data[x, y], this);
            image.Data[x, y] = new DColor(Colors[closest_color_dist_and_i.Item2]);
          }
        }
        //take double image data and save it back into the original Bitmap object
        image.SaveImageDataToBitmap();
      }
    };

    /// <summary>
    /// Generate a colormap based on input image.
    /// </summary>
    /// <param name="input">Input raster image.</param>
    /// <param name="desired_color_count">Required colormap size (ignore it if you must).</param>
    /// <param name="colors">Output palette (array of colors).</param>
    public static void Generate (Bitmap input, int desired_color_count, out Color[] colors)
    {
      DImage image = new DImage(input);
      //count how many colors are present in the image, if there are more than the palette can fit, set image_color_count to numCol + 1
      int found_color_count = image.CountColorsUntil(desired_color_count+1);

      //image has less or equal colors than the user wants in the palette
      bool image_too_few_colors = found_color_count <= desired_color_count;

      //how many colors will be in the palette - desired count (or less if image has less colors)
      int palette_color_count = image_too_few_colors ? found_color_count : desired_color_count;

      colors = new Color[palette_color_count];
      Palette best_palette = Palette.GetRandom(palette_color_count);

      //if there is less or equal amount of colors in the image than numCol, the palette are only the colors present in the image, get them
      if (image_too_few_colors)
      {
        best_palette.Colors = image.GetFirstNColors(palette_color_count);
      }
      //else, solve for palette the normal way
      else
      {
        ImageDistance image_dist = new ImageDistance(image);

        //these line were used to test accuracy of image distance methods
        //image_dist.DistanceApproxTestCube(image_color_count);
        //image_dist.DistanceApproxTest(image_color_count);

        double best_distance = double.PositiveInfinity;
        
        for (int i = 0; i < palette_test_count; i++)
        {
          Palette palette = new Palette(image, image_dist, palette_color_count);
          palette.SolveForPaletteColors(image, image_dist);
          double d = image_dist.DistanceFromImage(palette);
          Console.WriteLine(d);
          if (d < best_distance)
          {
            for (int j = 0; j < palette.ColorCount; j++)
            {
              best_palette.Colors[j] = new DColor(palette.Colors[j]);
            }
            best_distance = d;
          }
        }
        if (convert_image_to_palette)
        {
#pragma warning disable 162
          best_palette.ConvertImage(image, image_dist);
#pragma warning restore 162
        }
      }
      
      //convert computed double colors back to uint8 colors
      for (int i = 0; i < best_palette.ColorCount; i++)
      {
        colors[i] = best_palette.Colors[i].ToColor();
      }
    }
  }
}
