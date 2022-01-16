/**
 * Tries to simulate how rays bend around a black hole.
 */







// Default = contrast enhancing function.
formula.pixelTransform0 = (
  in ImageContext ic,
  ref float R,
  ref float G,
  ref float B) =>
{
  R = Util.Saturate(0.5f + 1.3f * (R - 0.5f));
  G = Util.Saturate(0.5f + 1.3f * (G - 0.5f));
  B = Util.Saturate(0.5f + 1.3f * (B - 0.5f));

  // Output color was modified.
  return true;
};


//class containing a 3D vector
class vec3
{
  double[] data;
  public double x
  {
    get => data[0];
    set => data[0] = value;
  }
  public double y
  {
    get => data[1];
    set => data[1] = value;
  }
  public double z
  {
    get => data[2];
    set => data[2] = value;
  }
  public vec3(double x, double y, double z)
  {
    data = new double[] { x, y, z };
  }
  public static vec3 operator *(vec3 v, double k)
  {
    return new vec3(v.x * k, v.y * k, v.z * k);
  }
  public double lengthSquared()
  {
    return x * x + y * y + z * z;
  }
  public double length()
  {
    return Math.Sqrt(lengthSquared());
  }
  public vec3 normalize()
  {
    double l = length();
    return new vec3(x / l, y / l, z / l);
  }
  public vec3 cross(vec3 v)
  {
    return new vec3(
        y * v.z - z * v.y,
        z * v.x - x * v.z,
        x * v.y - y * v.x);
  }
  public static vec3 operator +(vec3 a, vec3 b)
  {
    return new vec3(a.x + b.x, a.y + b.y, a.z + b.z);
  }
  public double dot(vec3 a)
  {
    return a.x * x + a.y * y + a.z * z;
  }
  public static vec3 operator *(vec3 a, vec3 b)
  {
    return new vec3(a.x * b.x, a.y * b.y, a.z * b.z);
  }
}




class Trace
{
  //all variables, although public, should only be modified during initialization
  public static double sphere_radius = 1;
  public static double light_speed = 2;
  public static double sphere_mass = 0.5;
  public static double time_step = 0.1;
  public static int max_iterations = 100;
  public static double ring_radius = 10;


  /**
   * Tracing - each step, move light particle according to current velocity. Then adjust velocity according to gravity of the black hole.
   * Check, whether the particle passed through an object(rings or the black sphere itself) during the last step
   * If yes, return interseciton point, otherwise continue with the next step
   */
  public static TracedT intersectAnyTraced(in vec3 ray_pos, in vec3 ray_dir)
  {
    //light particle position and velocity
    vec3 ray_p = ray_pos * 1.0;
    vec3 ray_v = ray_dir.normalize() * light_speed;
    for (int i = 0; i < max_iterations; i++)
    {
      vec3 ray = ray_v * time_step;
      RayT tg = intersectGround(ray_p, ray);
      RayT ts = intersectSphere(ray_p, ray);
      //if hit ground
      if (tg.isLimited())
      {
        //if during the last step, the particle hit both the ground and the sphere, return the one that was hit earlier
        if (ts.isLimited())
        {
          if (tg < ts)
          {
            return new TracedT(Hittable.GROUND, tg.ray(ray_p, ray));
          }
          else
          {
            return new TracedT(Hittable.SPHERE, ts.ray(ray_p, ray));
          }
        }
        //if only ground was hit, return it
        return new TracedT(Hittable.GROUND, tg.ray(ray_p, ray));
      }
      else
      {
        //if only sphere was hit, return it
        if (ts.isLimited())
        {
          return new TracedT(Hittable.SPHERE, ts.ray(ray_p, ray));
        }
      }
      //update particle position
      ray_p += ray_v * time_step;
      //update particle velocity - direction directly to the black hole (at (0, 0, 0)), proportional to time step and hole mass, inversely to the distance between
      ray_v += ray_p.normalize() * (-time_step * sphere_mass / ray_p.lengthSquared());
    }
    //return if nothing was hit
    return new TracedT(Hittable.NONE, new vec3(0, 0, 0));
  }

  //intersect given ray with ground
  public static RayT intersectGround(in vec3 ray_pos, in vec3 ray_dir)
  {
    //t of hit with ground
    double t = -ray_pos.z / ray_dir.z;
    //if ray travelled backwards or hit would happen too far away, return no hit
    if (t < 0 || Math.Abs(ray_dir.z) < 1e-10) return RayT.invalid;
    RayT x = new RayT(t);
    //if ground was hit outside of radius specified by rings, return no hit
    if (x.ray(ray_pos, ray_dir).length() > ring_radius) return RayT.invalid;
    //else return the hit
    return x;
  }

  //solving equations: x^2+y^2+z^2 = sphere_radius^2 && ray_pos + ray_dir * t = (x, y, z)
  //substituting into sphere equation, we can solve a quadratic equation for t:
  public static RayT intersectSphere(in vec3 ray_pos, in vec3 ray_dir)
  {
    //compute a, b, c coefficients for the quadratic equation
    double A = ray_dir.lengthSquared();
    double B = 2 * ray_pos.dot(ray_dir);
    double C = ray_pos.lengthSquared() - sphere_radius * sphere_radius;
    //compute discriminant
    double D = B * B - 4 * A * C;
    //if D<0, no hits
    if (D < 0) return RayT.invalid;
    //else compute t of the first hit and return it
    return new RayT((-B - Math.Sqrt(D)) / 2 / A);
  }
}


//holds t parameter in the '(x,y,z) = p + v * t' ray equation
class RayT
{
  double t;
  static double invalidValue = double.PositiveInfinity;

  public RayT(double t_) => t = t_;

  public static RayT invalid => new RayT(invalidValue);

  public bool isInvalid() => (t == invalidValue);

  public static bool operator <(RayT a, RayT b)
  {
    return a.t < b.t;
  }
  public static bool operator >(RayT a, RayT b)
  {
    return a.t > b.t;
  }
  public vec3 ray(in vec3 ray_pos, in vec3 ray_dir)
  {
    return ray_pos + ray_dir * t;
  }
  public bool isLimited()
  {
    return t >= 0 && t <= 1;
  }
}

//all hittable objects in the scene
public enum Hittable
{
  SPHERE,
  GROUND,
  NONE
}
//holds data about a hit done by a light particle
class TracedT
{
  public vec3 pos;
  public Hittable object_hit;
  public TracedT(Hittable obj, vec3 pos_)
  {
    object_hit = obj;
    pos = pos_;
  }
}

//these should be modified only during contextCreate!
class GlobalConsts
{
  public static vec3 cam_pos = new vec3(-10, -3, 1);
  public static vec3 cam_dir = new vec3(10, 0, -2).normalize();

  public static double cam_angle_modifier = 0.3;

  public static int sample_count = 5;
  public static (double x, double y)[] sample_positions;
}

//func for reading vec3 of given name from params
void tryReadVec3(Dictionary<string, string> dict, string name, ref vec3 target, vec3 default_val)
{
  if (dict.ContainsKey(name))
  {
    string[] values = dict[name].Split(' ');
    if (values.Length != 3)
    {
      Console.WriteLine("Invalid parameter value: {0}.", name);
      return;
    }
    try
    {
      target = new vec3(double.Parse(values[0]), double.Parse(values[1]), double.Parse(values[2]));
    }
    catch (System.FormatException)
    {
      Console.WriteLine("Invalid parameter format: {0}", name);
    }
  }
  else
  {
    target = default_val;
  }
}

//func for reading double of given name from params
void tryReadDouble(Dictionary<string, string> dict, string name, ref double target, double default_val)
{
  if (dict.ContainsKey(name))
  {
    try
    {
      target = double.Parse(dict[name]);
    }
    catch (System.FormatException)
    {
      Console.WriteLine("Invalid parameter format: {0}", name);
    }
  }
  else
  {
    target = default_val;
  }
}
//func for reading integer parameter of given name from params
void tryReadInt(Dictionary<string, string> dict, string name, ref int target, int default_val)
{
  if (dict.ContainsKey(name))
  {
    try
    {
      target = int.Parse(dict[name]);
    }
    catch (System.FormatException)
    {
      Console.WriteLine("Invalid parameter format: {0}", name);
    }
  }
  else
  {
    target = default_val;
  }
}

formula.contextCreate = (in Bitmap input, in string param) =>
{
  if (string.IsNullOrEmpty(param))
    return null;

  Dictionary<string, string> p = Util.ParseKeyValueList(param);


  //check whether any parameters are present
  tryReadVec3(p, "cam_pos", ref GlobalConsts.cam_pos, new vec3(-10, -3, 1));
  tryReadVec3(p, "cam_dir", ref GlobalConsts.cam_dir, new vec3(10, 0, -2));
  tryReadDouble(p, "cam_fov", ref GlobalConsts.cam_angle_modifier, 0.3);
  tryReadDouble(p, "sphere_r", ref Trace.sphere_radius, 1);
  tryReadDouble(p, "sphere_m", ref Trace.sphere_mass, 0.5);
  tryReadDouble(p, "light_speed", ref Trace.light_speed, 2);
  tryReadDouble(p, "time_step", ref Trace.time_step, 0.1);
  tryReadDouble(p, "ring_r", ref Trace.ring_radius, 10);
  tryReadInt(p, "trace_its", ref Trace.max_iterations, 100);
  tryReadInt(p, "sample_count", ref GlobalConsts.sample_count, 4);

  Random r = new Random();
  GlobalConsts.sample_positions = new (double, double)[GlobalConsts.sample_count];
  for (int i = 0; i < GlobalConsts.sample_count; i++)
  {
    GlobalConsts.sample_positions[i] = (r.NextDouble(), r.NextDouble());
  }

  //return tooltip info
  Dictionary<string, object> sc = new Dictionary<string, object>();

  sc["tooltip"] = "-------------------ANTIALIASING-------------------\n" +
                  "sample_count=<int> .. How many samples per pixel to use. Default is 4.\n" +
                  "-------------------BLACK HOLE RENDERER PARAMS-------------------\n" +
                  "cam_pos=<float> <float> <float> .. position of camera in world coordinates. Black hole is at (0, 0, 0). Formatting example:\"cam_pos=3 7 8\". Default = (-10, -3, 1)\r" +
                  "cam_dir=<float> <float> <float> .. direction of camera in world coordinates. Default = (10, 0, -2)\n" +
                  "cam_fov=<float> .. Between positive float, represents field of view. FOV = 90 * coeff degrees. Default 0.3(27 degrees)\n" +
                  "sphere_r=<float> .. sphere radius. Default = 1\n" +
                  "sphere_m=<float> .. sphere mass. Default = 0.5\n" +
                  "light_speed=<float> .. how fast the particles used to represent light will be moving. Default = 2\n" +
                  "time_step=<float> .. tracing time step. Default = 0.1\n" +
                  "ring_r=<double> .. radius of rings surrounding the black hole. Default = 10\n" +
                  "trace_its=<int> how many iterations will be used while tracing. Default = 100\n";

  return sc;
};

class RingsColor
{
  //mix computeColorCenter and computeColorCenter, create ring patterns
  public static vec3 computeColor(vec3 pos)
  {
    double d = pos.lengthSquared() / 25;
    //compute a simple sigmoid function (x / (1 + abs(x)), shifted by 1.5 to the right
    double din = d - 1.5;
    double w1 = (din / (1 + Math.Abs(din)) + 1) / 2;

    //use the weight value to mix the colors of the fringes and the center
    vec3 color = computeColorFringes(pos) * w1 + computeColorCenter(pos) * (1 - w1);

    //to make the pattern more ring-like, multiply the color with a sin function
    double k = Math.Pow(Math.Sin(Math.Pow(d, 1.2) * 10), 2);
    return color * k * new vec3(1, 0.5, 0.25);
  }

  //compute color at the center - this is mostly irregular, random mess of colors
  private static vec3 computeColorCenter(vec3 pos)
  {
    pos = pos * 0.2;

    double d = pos.lengthSquared();
    double ddiv = 1.0 / Math.Max(0.1, d - 0.5);

    //these functions dont really have any deep purpose
    //they are just chaotically nested sines, multiplications and additions to create mostly random noise

    pos.x += Math.Sin(pos.dot(new vec3(4, 4, 0))) * ddiv;
    pos.y += Math.Cos(pos.dot(new vec3(5, -7, 0)) + 2) * ddiv;
    pos.z += Math.Pow(Math.Sin(pos.dot(new vec3(-3, 0.5, 0)) + 8) * ddiv, 2);

    return new vec3(
        (Math.Sin(pos.dot(new vec3(5, 7, 2))) + 1),
        (Math.Sin(pos.dot(new vec3(2, -4, -6)) + 7) + 1),
        (Math.Sin(pos.dot(new vec3(8, -2, 4)) + 10) + 1)
    );
  }
  //compute color at the fringes - these are more continuous, but still chaotic, patterns
  private static vec3 computeColorFringes(vec3 pos)
  {
    pos = pos * 0.2;
    vec3 color_fringes = new vec3(
        Math.Sin(pos.dot(pos * new vec3(2, -5, 3)) * 5 + 1),
        Math.Sin(pos.dot(pos * new vec3(-3, -1, 7)) * 7 + 5),
        Math.Sin(pos.dot(pos * new vec3(8, 3, 5)) * 11 + 6)
    );
    color_fringes.x = Math.Sin(color_fringes.dot(new vec3(5, 4, 5.5)));
    color_fringes.y = Math.Sin(color_fringes.dot(new vec3(-2, -6, 2)));
    color_fringes.z = Math.Sin(color_fringes.dot(new vec3(1, 7, 2)));
    color_fringes.x = Math.Sin(color_fringes.y);
    return color_fringes;
  }
};

class BlackHole
{
  //directions up & right in camera viewport
  static vec3 cam_right = new vec3(0, 0, 1).cross(GlobalConsts.cam_dir).normalize() * GlobalConsts.cam_angle_modifier;
  static vec3 cam_up = cam_right.cross(GlobalConsts.cam_dir).normalize() * GlobalConsts.cam_angle_modifier;
  public static vec3 computeColor(double scr_x, double scr_y, double ratio)
  {
    //compute current ray direction
    vec3 ray = GlobalConsts.cam_dir.normalize() + cam_right * ratio * scr_x + cam_up * scr_y;

    //trace light particle to determine which object will be hit and where
    TracedT traced = Trace.intersectAnyTraced(GlobalConsts.cam_pos, ray);
    if (traced.object_hit == Hittable.NONE) //if nothing was hit, set background to very dark blue
    {
      return new vec3(0, 0, 0.1);
    }
    else if (traced.object_hit == Hittable.SPHERE)  //if sphere was hit, set background to pure black
    {
      return new vec3(0, 0, 0);
    }
    else      //if ground(=the rings) was hit, compute it's color 
    {
      return RingsColor.computeColor(traced.pos);
    }
  }
}

formula.pixelCreate = (
  in ImageContext ic,
  out float R,
  out float G,
  out float B) =>
{

  //ratio of image width to height - used when creating camera
  double ratio = 1.0 * ic.width / ic.height;

  // [x, y] in [-1, 1]


  vec3 color = new vec3(0, 0, 0);
  foreach ((double x, double y) sp in GlobalConsts.sample_positions)
  {
    double scr_x = (ic.x + sp.x) / (double)Math.Max(1, ic.width - 1) * 2 - 1;
    double scr_y = (ic.y + sp.y) / (double)Math.Max(1, ic.height - 1) * 2 - 1;
    color = color + BlackHole.computeColor(scr_x, scr_y, ratio);
  }
  color = color * (1 / (double)GlobalConsts.sample_count);

  R = (float)color.x;
  G = (float)color.y;
  B = (float)color.z;
};

