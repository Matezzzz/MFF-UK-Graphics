﻿using System;
using System.Drawing;
using System.Windows.Forms;

namespace _003colortransform
{
  public partial class FormColorTransform : Form
  {
    protected Bitmap inputImage  = null;
    protected Bitmap outputImage = null;

    static readonly string rev = "$Rev$".Split( ' ' )[ 1 ];

    public FormColorTransform ()
    {
      InitializeComponent();
      Text += " (rev: " + rev + ')';
    }

    private void buttonOpen_Click ( object sender, EventArgs e )
    {
      OpenFileDialog ofd = new OpenFileDialog();

      ofd.Title = "Open Image File";
      ofd.Filter = "Bitmap Files|*.bmp" +
          "|Gif Files|*.gif" +
          "|JPEG Files|*.jpg" +
          "|PNG Files|*.png" +
          "|TIFF Files|*.tif" +
          "|All image types|*.bmp;*.gif;*.jpg;*.png;*.tif";

      ofd.FilterIndex = 6;
      ofd.FileName = "";
      if ( ofd.ShowDialog() != DialogResult.OK )
        return;

      pictureBox1.Image = null;
      if ( inputImage != null )
        inputImage.Dispose();
      inputImage = (Bitmap)Image.FromFile( ofd.FileName );

      recompute();
    }

    private void recompute ()
    {
      if ( inputImage == null ) return;

      pictureBox1.Image = inputImage;
      if ( outputImage != null )
        outputImage.Dispose();

      Transform.TransformImage( inputImage, out outputImage, (double)numericParam.Value );

      pictureBox1.Image = outputImage;
    }

    private void buttonSave_Click ( object sender, EventArgs e )
    {
      if ( outputImage == null ) return;

      SaveFileDialog sfd = new SaveFileDialog();
      sfd.Title = "Save PNG file";
      sfd.Filter = "PNG Files|*.png";
      sfd.AddExtension = true;
      sfd.FileName = "";
      if ( sfd.ShowDialog() != DialogResult.OK )
        return;

      outputImage.Save( sfd.FileName, System.Drawing.Imaging.ImageFormat.Png );
    }

    private void buttonRedraw_Click ( object sender, EventArgs e )
    {
      recompute();
    }
  }
}
