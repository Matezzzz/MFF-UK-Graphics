using System;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.Threading;
using System.Windows.Forms;
using MathSupport;
using NCalc;
using OpenglSupport;
using OpenTK;
using OpenTK.Graphics.OpenGL;
using Utilities;

namespace _113graph
{
  public class Graph
  {
    /// <summary>
    /// Form-data initialization.
    /// </summary>
    public static void InitParams ( out string param, out string tooltip, out string expr, out string name, out MouseButtons trackballButton )
    {
      param = "domain=[-1.0;1.0;-1.0;1.0]";
      tooltip = "domain";
      expr = "1.0";
      trackballButton = MouseButtons.Left;

      name = "Matěj Mrázek";
    }

    /// <summary>
    /// Vertex array ([color] [normal] coord), index array
    /// </summary>
    uint[] VBOid = null;

    /// <summary>
    /// Currently allocated lengths of VBOs (in bytes)
    /// </summary>
    long[] VBOlen = null;

    /// <summary>
    /// Stride for vertex-array (in bytes).
    /// </summary>
    int stride = 0;

    /// <summary>
    /// Number of vertices (indices) to draw..
    /// </summary>
    int vertices = 0;

    /// <summary>
    /// Use vertex colors for rendering.
    /// </summary>
    bool useColors = true;

    /// <summary>
    /// Use normal vectors for rendering.
    /// </summary>
    bool useNormals = false;

    /// <summary>
    /// Cached expression.
    /// </summary>
    string expression = "";

    /// <summary>
    /// Cached param string.
    /// </summary>
    string param = "";

    /// <summary>
    /// Current model's center point.
    /// </summary>
    public Vector3 center = Vector3.Zero;

    /// <summary>
    /// Current model's diameter.
    /// </summary>
    public float diameter = 4.0f;

    /// <summary>
    /// Near point for the current model.
    /// </summary>
    public float near = 0.1f;

    /// <summary>
    /// Far point for the current model.
    /// </summary>
    public float far = 20.0f;

    public void InitOpenGL ( GLControl glc )
    {
      // log OpenGL info just for curiosity:
      GlInfo.LogGLProperties();

      // general OpenGL:
      glc.VSync = true;
      GL.ClearColor( Color.FromArgb( 14, 20, 40 ) );    // darker "navy blue"
      GL.Enable( EnableCap.DepthTest );
      GL.Enable( EnableCap.VertexProgramPointSize );
      //GL.Enable()
      GL.ShadeModel( ShadingModel.Smooth );

      // VBO init:
      VBOid = new uint[ 2 ];           // one big buffer for vertex data, another buffer for tri/line indices
      GL.GenBuffers( 2, VBOid );
      GlInfo.LogError( "VBO init" );
      VBOlen = new long[ 2 ];
    }

    public void InitSimulation ( string par, string expr )
    {
      param = "";
      expression = "";
      RegenerateGraph( par, expr );
    }

    /// <summary>
    /// Recompute the graph, prepare VBO content and upload it to the GPU...
    /// </summary>
    /// <param name="param">New param string.</param>
    /// <param name="expr">New expression string.</param>
    /// <returns>null if OK, error message otherwise.</returns>
    public string RegenerateGraph ( string par, string expr )
    {
      // !!!{{ TODO: add graph data regeneration code here

      if ( expr == expression &&
           par == param )
        return null;                // nothing to do..

      double xMin = -1.0, xMax = 1.0, yMin = -1.0, yMax = 1.0;

      // input params:
      Dictionary<string, string> p = Util.ParseKeyValueList( param );

      // domain: [xMin;xMax;yMin;yMax]
      List<double> dom = null;
      if ( Util.TryParse( p, "domain", ref dom, ';' ) &&
           dom != null &&
           dom.Count >= 4 )
      {
        xMin = dom[ 0 ];
        xMax = Math.Max( xMin + 1.0e-6, dom[ 1 ] );
        yMin = dom[ 2 ];
        yMax = Math.Max( yMin + 1.0e-6, dom[ 3 ] );
      }

      int resolution_x = 100;
      int resolution_y = 100;

      // Data for VBO:
      useColors = true;
      useNormals = false;
      stride = Vector3.SizeInBytes * (1 + (useColors ? 1 : 0) + (useNormals ? 1 : 0));
      long newVboSize = stride * (resolution_x + 1) * (resolution_y + 1);     // pilot .. three vertices
      vertices = 6 * resolution_x * resolution_y;                     // pilot .. three indices
      long newIndexSize = sizeof( uint ) * vertices;

      // OpenGL stuff
      GL.EnableClientState( ArrayCap.VertexArray );
      if ( useColors )
        GL.EnableClientState( ArrayCap.ColorArray );
      if ( useNormals )
        GL.EnableClientState( ArrayCap.NormalArray );

      // Vertex array: [color] [normal] coordinate
      GL.BindBuffer( BufferTarget.ArrayBuffer, VBOid[ 0 ] );
      if ( newVboSize != VBOlen[ 0 ] )
      {
        VBOlen[ 0 ] = newVboSize;
        GL.BufferData( BufferTarget.ArrayBuffer, (IntPtr)VBOlen[ 0 ], IntPtr.Zero, BufferUsageHint.DynamicDraw );
      }

      IntPtr videoMemoryPtr = GL.MapBuffer( BufferTarget.ArrayBuffer, BufferAccess.WriteOnly );


      double func_min = double.PositiveInfinity;
      double func_max = double.NegativeInfinity;

      Expression e = new Expression(expr);

      double scale(int i, int steps, double min, double max) => min + (max - min) * i / steps;

      double[,,] function_values = new double[resolution_x+1,resolution_y+1,3];
      for (int x_i = 0; x_i <= resolution_x; x_i++)
      {
        double x = scale(x_i, resolution_x, xMin, xMax);
        for (int y_i = 0; y_i <= resolution_y; y_i++)
        {
          double y = scale(y_i, resolution_y, yMin, yMax);
          
          double result;
          try
          {
            e.Parameters["x"] = x;
            e.Parameters["y"] = y;
            e.Parameters["z"] = y;
            result = (double)e.Evaluate();
            function_values[x_i, y_i, 0] = x;
            function_values[x_i, y_i, 1] = result;
            function_values[x_i, y_i, 2] = y;

            func_min = Math.Min(func_min, result);
            func_max = Math.Max(func_max, result);

            if (double.IsNaN(result))
              throw new Exception("NCalc: NaN");
          }
          catch (Exception ex)
          {
            return ex.Message;
          }
        }
      }




      double ReLU (double x) => (x < 0) ? 0 : x;

      unsafe
      {
        float* ptr = (float*)videoMemoryPtr.ToPointer();

        for (int x_i = 0; x_i <= resolution_x; x_i++)
        {
          for (int y_i = 0; y_i <= resolution_y; y_i++)
          {
            double y_norm = (function_values[x_i, y_i, 1] - func_min) / (func_max - func_min);

            *ptr++ = (float)(1 - ReLU(2 * y_norm - 1));
            *ptr++ = (float)(1 - ReLU(1 - 2 * y_norm));
            *ptr++ = (float)0;

            *ptr++ = (float)function_values[x_i, y_i, 0];
            *ptr++ = (float)function_values[x_i, y_i, 1];
            *ptr++ = (float)function_values[x_i, y_i, 2];
          }
        }
      }
      
      // Everything seems to be OK:
      expression = expr;
      param = par;

      Vector3  v = new Vector3( (float) scale(1, 2, xMin, xMax), (float)scale(1, 2, yMin, yMax), (float)scale(1, 2, yMin, yMax));

      GL.UnmapBuffer( BufferTarget.ArrayBuffer );
      GL.BindBuffer( BufferTarget.ArrayBuffer, 0 );

      // Index buffer
      GL.BindBuffer( BufferTarget.ElementArrayBuffer, VBOid[ 1 ] );
      if ( newIndexSize != VBOlen[ 1 ] )
      {
        VBOlen[ 1 ] = newIndexSize;
        GL.BufferData( BufferTarget.ElementArrayBuffer, (IntPtr)VBOlen[ 1 ], IntPtr.Zero, BufferUsageHint.StaticDraw );
      }

      videoMemoryPtr = GL.MapBuffer( BufferTarget.ElementArrayBuffer, BufferAccess.WriteOnly );
      unsafe
      {
        uint coord(uint x, uint y)
        {
          return (uint) (y * (resolution_x + 1) + x);
        }
        uint* ptr = (uint*)videoMemoryPtr.ToPointer();
        for (uint x = 0; x < resolution_x; x++)
        {
          for (uint y = 0; y < resolution_y; y++)
          {
            *ptr++ = coord(x, y);
            *ptr++ = coord(x+1, y);
            *ptr++ = coord(x+1, y+1);

            *ptr++ = coord(x, y);
            *ptr++ = coord(x+1, y+1);
            *ptr++ = coord(x, y+1);
          }
        }
      }
      GL.UnmapBuffer( BufferTarget.ElementArrayBuffer );
      GL.BindBuffer( BufferTarget.ElementArrayBuffer, 0 );

      // Change the graph dimension:
      center   =     v;
      diameter =  2.0f;
      near     =  0.1f;
      far      = 20.0f;

      Form1.form.SetStatus( $"Tri: {vertices / 3}" );

      return null;
    }

    /// <summary>
    /// Rendering code itself (separated for clarity).
    /// </summary>
    public void RenderScene ( ref long primitiveCounter )
    {
      // Scene rendering:
      if ( Form1.form.drawGraph &&
           VBOlen[ 0 ] > 0L )        // buffers are nonempty => render
      {
        // [color] [normal] coordinate
        GL.BindBuffer( BufferTarget.ArrayBuffer, VBOid[ 0 ] );
        IntPtr p = IntPtr.Zero;
        if ( useColors )             // are colors present?
        {
          GL.ColorPointer( 3, ColorPointerType.Float, stride, p );
          p += Vector3.SizeInBytes;
        }
        if ( useNormals )            // are normals present?
        {
          GL.NormalPointer( NormalPointerType.Float, stride, p );
          p += Vector3.SizeInBytes;
        }
        GL.VertexPointer( 3, VertexPointerType.Float, stride, p );

        // index buffer
        GL.BindBuffer( BufferTarget.ElementArrayBuffer, VBOid[ 1 ] );

        // Multiple instancing of the scene:
        GL.DrawElements( PrimitiveType.Triangles, vertices, DrawElementsType.UnsignedInt, IntPtr.Zero );

        primitiveCounter += vertices / 3;
      }
      else                           // color cube
      {
        GL.Begin( PrimitiveType.Quads );

        GL.Color3( 0.0f, 1.0f, 0.0f );          // Set The Color To Green
        GL.Vertex3( 1.0f, 1.0f, -1.0f );        // Top Right Of The Quad (Top)
        GL.Vertex3( -1.0f, 1.0f, -1.0f );       // Top Left Of The Quad (Top)
        GL.Vertex3( -1.0f, 1.0f, 1.0f );        // Bottom Left Of The Quad (Top)
        GL.Vertex3( 1.0f, 1.0f, 1.0f );         // Bottom Right Of The Quad (Top)

        GL.Color3( 1.0f, 0.5f, 0.0f );          // Set The Color To Orange
        GL.Vertex3( 1.0f, -1.0f, 1.0f );        // Top Right Of The Quad (Bottom)
        GL.Vertex3( -1.0f, -1.0f, 1.0f );       // Top Left Of The Quad (Bottom)
        GL.Vertex3( -1.0f, -1.0f, -1.0f );      // Bottom Left Of The Quad (Bottom)
        GL.Vertex3( 1.0f, -1.0f, -1.0f );       // Bottom Right Of The Quad (Bottom)

        GL.Color3( 1.0f, 0.0f, 0.0f );          // Set The Color To Red
        GL.Vertex3( 1.0f, 1.0f, 1.0f );         // Top Right Of The Quad (Front)
        GL.Vertex3( -1.0f, 1.0f, 1.0f );        // Top Left Of The Quad (Front)
        GL.Vertex3( -1.0f, -1.0f, 1.0f );       // Bottom Left Of The Quad (Front)
        GL.Vertex3( 1.0f, -1.0f, 1.0f );        // Bottom Right Of The Quad (Front)

        GL.Color3( 1.0f, 1.0f, 0.0f );          // Set The Color To Yellow
        GL.Vertex3( 1.0f, -1.0f, -1.0f );       // Bottom Left Of The Quad (Back)
        GL.Vertex3( -1.0f, -1.0f, -1.0f );      // Bottom Right Of The Quad (Back)
        GL.Vertex3( -1.0f, 1.0f, -1.0f );       // Top Right Of The Quad (Back)
        GL.Vertex3( 1.0f, 1.0f, -1.0f );        // Top Left Of The Quad (Back)

        GL.Color3( 0.0f, 0.0f, 1.0f );          // Set The Color To Blue
        GL.Vertex3( -1.0f, 1.0f, 1.0f );        // Top Right Of The Quad (Left)
        GL.Vertex3( -1.0f, 1.0f, -1.0f );       // Top Left Of The Quad (Left)
        GL.Vertex3( -1.0f, -1.0f, -1.0f );      // Bottom Left Of The Quad (Left)
        GL.Vertex3( -1.0f, -1.0f, 1.0f );       // Bottom Right Of The Quad (Left)

        GL.Color3( 1.0f, 0.0f, 1.0f );          // Set The Color To Violet
        GL.Vertex3( 1.0f, 1.0f, -1.0f );        // Top Right Of The Quad (Right)
        GL.Vertex3( 1.0f, 1.0f, 1.0f );         // Top Left Of The Quad (Right)
        GL.Vertex3( 1.0f, -1.0f, 1.0f );        // Bottom Left Of The Quad (Right)
        GL.Vertex3( 1.0f, -1.0f, -1.0f );       // Bottom Right Of The Quad (Right)

        GL.End();

        primitiveCounter += 12;
      }
    }

    public void Intersect ( ref Vector3d p0, ref Vector3d p1, ref double nearest )
    {
      // Compute intersection of the given ray (p0, p1) with the scene,
      // act upon that (i.e. set status string, ..)
      Vector2d uv;

#if false

      Vector3 A, B, C;
      double curr = Geometry.RayTriangleIntersection( ref p0, ref p1, ref A, ref B, ref C, out uv );
      if ( !double.IsInfinity( curr ) &&
           curr < nearest )
        nearest = curr;

#else

      Vector3d ul   = new Vector3d( -1.0, -1.0, -1.0 );
      Vector3d size = new Vector3d( 2.0, 2.0, 2.0 );
      if ( Geometry.RayBoxIntersection( ref p0, ref p1, ref ul, ref size, out uv ) )
      {
        nearest = uv.X;
        Form1.form.SetStatus( string.Format( CultureInfo.InvariantCulture, "[{0:f2},{1:f2},{2:f2}]",
                                             p0.X + nearest * p1.X,
                                             p0.Y + nearest * p1.Y,
                                             p0.Z + nearest * p1.Z ) );
      }

#endif
    }

    /// <summary>
    /// Handles mouse-button push.
    /// </summary>
    /// <returns>True if handled.</returns>
    public bool MouseButtonDown ( MouseEventArgs e )
    {
      return false;
    }

    /// <summary>
    /// Handles mouse-button release.
    /// </summary>
    /// <returns>True if handled.</returns>
    public bool MouseButtonUp ( MouseEventArgs e )
    {
      return false;
    }

    /// <summary>
    /// Handles mouse move.
    /// </summary>
    /// <returns>True if handled.</returns>
    public bool MousePointerMove ( MouseEventArgs e )
    {
      return false;
    }

    /// <summary>
    /// Handles keyboard key press.
    /// </summary>
    /// <returns>True if handled.</returns>
    public bool KeyHandle ( KeyEventArgs e )
    {
      return false;
    }

    public void Destroy ()
    {
      if ( VBOid != null &&
           VBOid[ 0 ] != 0 )
      {
        GL.DeleteBuffers( 2, VBOid );
        VBOid = null;
      }
    }
  }
}
