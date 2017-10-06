using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;
using System.Globalization;

namespace Drawing
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

       
        static int n, m, h;

        Bitmap imageE, imageB;
        double[,] masE;
        double[,] masB;


        private void button1_Click(object sender, EventArgs e)
        {



            Stream myStream = null;
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = "\\...";
            openFileDialog1.Filter = "txt files (*.txt)|*.txt";
            openFileDialog1.FilterIndex = 2;
            openFileDialog1.RestoreDirectory = true;

            if (openFileDialog1.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    if ((myStream = openFileDialog1.OpenFile()) != null)
                    {
                        using (myStream)
                        {
                            StreamReader file = File.OpenText(@openFileDialog1.FileName);

                            n = Int32.Parse(file.ReadLine());
                            m = Int32.Parse(file.ReadLine());
                            h = Int32.Parse(file.ReadLine());
                            masE = new double[n, m];
                            masB = new double[n, m];
                            string[] str = new string[n * m];

                            for (int i = 0; i < n; i++)
                                for (int j = 0; j < m; j++)
                                {
                                    masE[i, j] = Double.Parse(file.ReadLine());
                                    masB[i, j] = Double.Parse(file.ReadLine());
                                }
                        }
                    }
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Error: Could not read file from disk. Original error: " + ex.Message);
                }
            }

            imageB = new Bitmap(n, m);
            imageE = new Bitmap(n, m);


            //StreamReader file = File.OpenText(@"input.txt");
            //for (int i = 0; i < n; i++)
            //    for (int j = 0; j < m; j++)
            //    {
            //        mas[i, j] = Double.Parse(file.ReadLine());
            //    }

            //using (System.IO.StreamWriter output =
            //new System.IO.StreamWriter(@"output.txt"))
            //{
            //    for (int i = 0; i < n; i++)
            //        for (int j = 0; j < m; j++)
            //        {
            //            output.WriteLine(mas[i, j]);
            //        }

            //}
            double min = masE[0, 0];
            double max = masE[0, 0];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                {
                    if (masE[i, j] > max)
                        max = masE[i, j];
                    if (masE[i, j] < min)
                        min = masE[i, j];
                }

            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                {
                    int col = (int)((masE[i, j] - min) / (max - min) * 255);
                   imageE.SetPixel(i, j, Color.FromArgb(col, col, col));
                }

            pictureBox1.Image = imageE;
            pictureBox1.Refresh();

            min = masB[0, 0];
            max = masB[0, 0];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                {
                    if (masB[i, j] > max)
                        max = masB[i, j];
                    if (masB[i, j] < min)
                        min = masB[i, j];
                }

            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                {
                    int col = (int)((masB[i, j] - min) / (max - min) * 255);
                    imageB.SetPixel(i, j, Color.FromArgb(col, col, col));
                }

            pictureBox2.Image = imageB;
            pictureBox2.Refresh();
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }

        private void pictureBox1_Click(object sender, EventArgs e)
        {

        }
    }
}
