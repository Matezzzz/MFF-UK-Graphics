﻿using System;
using System.Drawing;
using System.Windows.Forms;

namespace _006warping
{
  public partial class Form1 : Form
  {
    protected Bitmap inputImage  = null;
    protected Bitmap outputImage = null;

    public Form1 ()
    {
      InitializeComponent();
      String[] tok = "$Rev$".Split( ' ' );
      Text += " (rev: " + tok[ 1 ] + ')';
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

      pictureResult.SetPicture( null );
      if ( inputImage != null )
        inputImage.Dispose();
      inputImage = (Bitmap)Image.FromFile( ofd.FileName );

      recompute();
    }

    private void recompute ()
    {
      if ( inputImage == null ) return;

      pictureResult.SetPicture( inputImage );
    }

    private void buttonSave_Click ( object sender, EventArgs e )
    {
      if ( inputImage == null ) return;

      SaveFileDialog sfd = new SaveFileDialog();
      sfd.Title = "Save PNG file";
      sfd.Filter = "PNG Files|*.png";
      sfd.AddExtension = true;
      sfd.FileName = "";
      if ( sfd.ShowDialog() != DialogResult.OK )
        return;

      pictureResult.GetPicture().Save( sfd.FileName, System.Drawing.Imaging.ImageFormat.Png );
    }

    private void numericParam_ValueChanged ( object sender, EventArgs e )
    {
      recompute();
    }
  }
}
