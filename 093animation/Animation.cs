using System;
using System.Collections.Generic;
using System.Drawing;
using LineCanvas;
using MathSupport;
using Utilities;

namespace _093animation
{

  public static class ColorExtensions
  {
    public static Color Mult(this Color c, double k)
    {
      return Color.FromArgb((int) (k * c.R), (int)(k * c.G), (int)(k * c.B));
    }
    public static Color Add(this Color a, Color b)
    {
      return Color.FromArgb(a.R + b.R, a.G + b.G, a.B + b.B);
    }
  }

  class GloriousDraw
  {
    static double percent (double step, double step_count) => step / step_count;
    static double angle (double step, double step_count, double shift = 0) => 2 * Math.PI * percent(step + shift, step_count);
    static double color (double step, double step_count, double shift) => 360 * percent(step + shift, step_count);


    static void drawLine (Canvas c, double x, double y, double a1, double a2, double r1, double r2)
    {
      c.Line(x + Math.Sin(a1) * r1, y + Math.Cos(a1) * r1, x + Math.Sin(a2) * r2, y + Math.Cos(a2) * r2);
    }
    public static void DrawCircle (Canvas c, double x, double y, int steps, double r1, int angle_shift_steps, double start_angle, Color color)
    {
      for (int k = 0; k < steps; k++)
      {
        c.SetColor(color);
        double a = angle(k, steps, start_angle);
        drawLine(c, x, y, a, angle(k, steps, angle_shift_steps+start_angle), r1, r1);
        drawLine(c, x, y, a, angle(k, steps, start_angle+1), r1, r1);
      }
    }

    public static void DrawSpiral (Canvas c, double x, double y, int steps, double rdif, int color_hue_shift, int angle_shift_steps, int init_angle_steps, double color_lightness)
    {
      for (int k = 0; k < 10 * steps; k++)
      {
        c.SetColor(Arith.HSVToColor(color(k, steps, color_hue_shift), 1.0, color_lightness));
        int ang_k = init_angle_steps + k;
        double a = angle(ang_k, steps);
        drawLine(c, x, y, a, angle(ang_k, steps, angle_shift_steps), k + rdif, k);
        drawLine(c, x, y, a, angle(ang_k, steps, -angle_shift_steps), k + rdif, k);
      }
    }
  }

  public class Animation
  {
    /// <summary>
    /// Form data initialization.
    /// </summary>
    /// <param name="name">Your first-name and last-name.</param>
    /// <param name="wid">Image width in pixels.</param>
    /// <param name="hei">Image height in pixels.</param>
    /// <param name="from">Animation start in seconds.</param>
    /// <param name="to">Animation end in seconds.</param>
    /// <param name="fps">Frames-per-seconds.</param>
    /// <param name="param">Optional text to initialize the form's text-field.</param>
    /// <param name="tooltip">Optional tooltip = param help.</param>
    public static void InitParams (out string name, out int wid, out int hei, out double from, out double to, out double fps, out string param, out string tooltip)
    {
      // {{

      // Put your name here.
      name = "Matěj Mrázek";

      // Image size in pixels.
      wid = 1920;
      hei = 1080;

      // Animation.
      from = 0.0;
      to   = 20.0;
      fps  = 25.0;

      // Specific animation params.
      param = "anti=true,objects=10,trail_length=50,trail_resolution=3, gravity_coefficient=15000000, sphere_count=20, border_force=10";


      // Tooltip = help.
      tooltip = "anti=<int, 1 or 0>, objects=<int>, trail_length=<int>, trail_resolution=<int>";

      // }}
    }

    class Sphere
    {
      public double mass;
      public double radius;
      Color color;
      Color tail_color;

      double angle;
      double vel_angle;

      public double pos_x, pos_y;
      public double vel_x, vel_y;
      double force_x, force_y;


      public static int tail_length = 50;
      double tail_width;

      List<(double x, double y, double angle, Color color)> history;

      public Sphere (double pos_x_, double pos_y_, double mass_, double radius_, int r, int g, int b)
      {
        mass = mass_;
        radius = radius_;
        angle = 0;
        vel_angle = 1;
        pos_x = pos_x_;
        pos_y = pos_y_;
        vel_x = 0;
        vel_y = 0;
        force_x = 0;
        force_y = 0;

        tail_width = radius;

        color = Color.FromArgb(r, g, b);
        tail_color = Color.FromArgb(r / 2, g / 2, b / 2);

        history = new List<(double, double, double, Color)>();
      }
      public Sphere SetVelocity(double x, double y)
      {
        vel_x = x;
        vel_y = y;
        return this;
      }
      public void Move(double x, double y)
      {
        pos_x += x;
        pos_y += y;
      }

      public void AddForce(double x, double y)
      {
        force_x += x; force_y += y;
      }
      public void MixColors(Sphere b)
      {
        Color col1 = color;
        Color col2 = b.color;
        color   = col1.Mult(0.9).Add(col2.Mult(0.1));
        b.color = col2.Mult(0.9).Add(col1.Mult(0.1));
      }
      public void Update (double dt)
      {
        vel_x += force_x * dt / mass;
        vel_y += force_y * dt / mass;

        pos_x += vel_x * dt;
        pos_y += vel_y * dt;

        vel_x *= 0.999;
        vel_y *= 0.999;

        force_x = 0;
        force_y = 0;

        vel_angle = (Math.Abs(vel_x) + Math.Abs(vel_y)) / 50 + 2;
        angle += vel_angle * dt;
      }
      public void Save ()
      {
        history.Add((pos_x, pos_y, angle, color));
      }

      public void Draw (Canvas c, int frame_i)
      {
        c.SetColor(tail_color);

        bool drawing_enabled = false;


        int i_start = frame_i - tail_length;
        if (i_start < 1)
          i_start = 1;


        (double x, double y, double ang, Color col) = history[frame_i];
        (double prev_x1, double prev_y1, double _, Color _) = history[frame_i];
        (double prev_x2, double prev_y2, double _, Color _) = history[frame_i];
        for (int i = frame_i; i >= i_start; i--)
        {
          (double prevx, double prevy, double _, Color _) = history[i-1];
          (double curx, double cury, double _, Color _) = history[i];
          (double nextx, double nexty, double _, Color _) = history[i+1];

          double dist = Math.Sqrt(Math.Pow(curx - x, 2) + Math.Pow(cury - y, 2));
          if (dist > radius)
            drawing_enabled = true;

          double len_part = (1.0 + i - frame_i + tail_length) / tail_length;

          double dx = nextx - prevx; double dy = nexty - prevy;
          double lnx = -dy; double lny = dx;

          double ln = Math.Sqrt(lnx * lnx + lny * lny);

          lnx *= len_part / ln * tail_width;
          lny *= len_part / ln * tail_width;

          double x1 = curx + lnx, y1 = cury + lny, x2 = curx - lnx, y2 = cury - lny;


          if (drawing_enabled)
          {
            c.Line(prev_x1, prev_y1, x2, y2);
            c.Line(prev_x2, prev_y2, x1, y1);
          }

          prev_x1 = x1;
          prev_x2 = x2;
          prev_y1 = y1;
          prev_y2 = y2;
        }

        

        GloriousDraw.DrawCircle(c, x, y, (int) radius, radius, (int) (0.3*radius), ang, col);
       
      }
    }

    class GravitySim
    {
      Sphere[] spheres;


      public static double gravity_coef = 15000000;
      public static double border_force_k = 10;
      public static int sphere_count = 10;
      public static int trail_resolution = 3;

      public GravitySim (int wid, int hei)
      {
        spheres = new Sphere[sphere_count];

        spheres[0] = new Sphere(wid / 2.0, hei / 2.0, 1500, 50, 255, 0, 0);
        //spheres[1] = new Sphere(wid * 0.5, hei * 0.75, 10, 20, 0, 255, 0);

        
        Random rand = new Random();

        for (int i = 1; i < sphere_count; i++)
        {
          double a = rand.NextDouble() * 2 * Math.PI;
          double r = rand.NextDouble() * wid / 8 + wid / 4;
          double v = r;
          double x = Math.Sin(a), y = Math.Cos(a);
          double m = rand.Next(5, 25);
          spheres[i] = new Sphere(wid / 2 + x * r, hei / 2 + y * r, Math.Sqrt(m) * 2, m, rand.Next(255), rand.Next(255), rand.Next(255)).SetVelocity(-y * v, x * v);
          //double m = rand.Next(5, 25);
          //spheres[i] = new Sphere(rand.Next(wid), rand.Next(hei), m, Math.Sqrt(m) * 4, rand.Next(255), rand.Next(255), rand.Next(255)).SetVelocity(rand.Next(-200, 200), rand.Next(-200, 200));
        }

      }
      public void Run (int frame_count, double dt_frame, int wid, int hei)
      {
        double dt = dt_frame / trail_resolution;
        for (int f = 0; f < frame_count * trail_resolution; f++)
        {
          for (int i = 0; i < spheres.Length; i++)
          {
            Sphere a = spheres[i];
            for (int j = i + 1; j < spheres.Length; j++)
            {
              Sphere b = spheres[j];
              double dif_x = a.pos_x - b.pos_x;
              double dif_y = a.pos_y - b.pos_y;
              double dist = Math.Sqrt(dif_x * dif_x + dif_y * dif_y);

              double M = gravity_coef * a.mass * b.mass * dt / (dist * dist * dist);

              //prevent glitches caused by forces being too large
              const int M_max = 1000;
              if (M > M_max)
                M = M_max;
              double f_x = M * dif_x;
              double f_y = M * dif_y;
              a.AddForce(-f_x, -f_y);
              b.AddForce(f_x, f_y);

              double intersection = a.radius + b.radius - dist;
              if (intersection > 0)
              {
                double d2 = intersection / (b.mass / a.mass + 1);
                double d1 = intersection - d2;
                d2 /= dist;
                d1 /= dist;

                a.Move(d1 * dif_x, d1 * dif_y);
                b.Move(-d2 * dif_x, -d2 * dif_y);


                double crash_en = ((a.vel_x * a.mass - b.vel_x * b.mass) * dif_x + (a.mass * a.vel_y - b.mass * b.vel_y) * dif_y) / (dist * dist);
                double j_x = crash_en * dif_x / dt;
                double j_y = crash_en * dif_y / dt;
                a.AddForce(-j_x, -j_y);
                b.AddForce(j_x, j_y);
                a.MixColors(b);
              }
            }
            
            if (a.pos_x < 0)
              a.AddForce(-border_force_k * a.pos_x, 0);
            if (a.pos_x > wid)
              a.AddForce(-border_force_k * (a.pos_x - wid), 0);
            if (a.pos_y < 0)
              a.AddForce(0, -border_force_k * a.pos_y);
            if (a.pos_y > hei)
              a.AddForce(0, -border_force_k * (a.pos_y - hei));
          }
          foreach (Sphere s in spheres)
          {
            s.Save();
            s.Update(dt);
          }
        }
      }
      public void Draw(Canvas c, int frame)
      {
        foreach (Sphere s in spheres)
        {
          s.Draw(c, frame * trail_resolution);
        }
      }
    }

    /// <summary>
    /// Global initialization. Called before each animation batch
    /// or single-frame computation.
    /// </summary>
    /// <param name="width">Width of the future canvas in pixels.</param>
    /// <param name="height">Height of the future canvas in pixels.</param>
    /// <param name="start">Start time (t0)</param>
    /// <param name="end">End time (for animation length normalization).</param>
    /// <param name="fps">Required fps.</param>
    /// <param name="param">Optional string parameter from the form.</param>
    ///

    static GravitySim simulation;
    static double fps;

    public static void InitAnimation (int width, int height, double start, double end, double fps_, string param)
    {
      Dictionary<string, string> user_params = Util.ParseKeyValueList(param);
      Util.TryParseIntRange(user_params, "trail_resolution", ref GravitySim.trail_resolution);
      Util.TryParseRational(user_params, "gravity_coefficient", ref GravitySim.gravity_coef);
      Util.TryParseIntRange(user_params, "sphere_count", ref GravitySim.sphere_count);
      Util.TryParseRational(user_params, "border_force", ref GravitySim.border_force_k);
      Util.TryParseIntRange(user_params, "trail_length", ref Sphere.tail_length);


      fps = fps_;
      simulation = new GravitySim(width, height);
      simulation.Run((int)((end - start) * fps) + 100, 1.0 / fps, width, height);

      
  }

    /// <summary>
    /// Draw single animation frame.
    /// Has to be re-entrant!
    /// </summary>
    /// <param name="c">Canvas to draw to.</param>
    /// <param name="time">Current time in seconds.</param>
    /// <param name="start">Start time (t0)</param>
    /// <param name="end">End time (for animation length normalization).</param>
    /// <param name="param">Optional string parameter from the form.</param>
    public static void DrawFrame (Canvas c, double time, double start, double end, string param)
    {
      Dictionary<string, string> user_params = Util.ParseKeyValueList(param);
      int anti = 1;
      Util.TryParseIntRange(user_params, "anti", ref anti);
      c.SetAntiAlias(anti != 0);
      simulation.Draw(c, (int) (time * fps));
    }
  }
}
