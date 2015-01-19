using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Diagnostics;
using System.Drawing;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using RTree;
using Point = RTree.Point;
using Rectangle = RTree.Rectangle;

namespace mbR_NGQ
{
    public partial class RTreeViewer : Form
    {
        private RTree<Point> rTree;

        public RTreeViewer()
        {
            InitializeComponent();
            
            System.Diagnostics.Debug.WriteLine("Data Generation Start");

            DataSet dataSet = new DataSet();
            dataSet.DataGenerator = new UniformGenerator();
            dataSet.DataGeneration(Config.n * Config.m);

            System.Diagnostics.Debug.WriteLine("Data Generation End");
            System.Diagnostics.Debug.WriteLine("RTree Insert Start");

            rTree = new RTree<Point>();
            foreach (Point p in dataSet.Points)
            {
                rTree.Add(new Rectangle(p.coordinates, p.coordinates), p);
            }
            rTree.CalculateBitArray();

            System.Diagnostics.Debug.WriteLine("RTree Insert End");

            List<RTree<Point>.MGroup<Point>> nearGroups = null;

            Stopwatch sw = new Stopwatch();
            sw.Start();

            Random rand = new Random();
            Point queryPoint = null;

            for(int i = 0; i < 100; i++)
            {
                float x = (float)(rand.NextDouble() * (Config.maxX - Config.minX) + Config.minX);
                float y = (float)(rand.NextDouble() * (Config.maxY - Config.minY) + Config.minY);
                queryPoint = new Point((float)x, (float)y, 0);
                nearGroups = rTree.NearestGroup(queryPoint);
                //System.Diagnostics.Debug.Assert(rTree.CheckNearGroup(queryPoint, nearGroups) == true);
            }
            
            /*float x = (float)(rand.NextDouble() * (Config.maxX - Config.minX) + Config.minX);
            float y = (float)(rand.NextDouble() * (Config.maxY - Config.minY) + Config.minY);
            queryPoint = new Point((float)x, (float)y, 0);
            nearGroups = rTree.NearestGroup(queryPoint);*/

            sw.Stop();
            System.Diagnostics.Debug.WriteLine(sw.ElapsedMilliseconds);

            /*rTreePanel1.RTree = rTree;*/
            rTreePanel1.NearGroups = nearGroups;
            rTreePanel1.QueryPoint = queryPoint;
            rTreePanel1.InitRTree();

            SetLevelList();
            SetGroupList();
        }

        private void SetLevelList()
        {
            this.checkedListBox1.Items.Clear();
            for (int i = 0; i < rTree.TreeHeight; i++)
            {
                this.checkedListBox1.Items.Add(string.Format("level {0}", i+1));
            }
        }

        private void SetGroupList()
        {
            this.checkedListBox2.Items.Clear();
            for (int i = 0; i < rTreePanel1.NearGroups.Count; i++)
            {
                this.checkedListBox2.Items.Add(string.Format("Group {0}", i + 1));
            }
        }

        private void checkedListBox1_ItemCheck(object sender, ItemCheckEventArgs e)
        {
            if(e.NewValue == CheckState.Checked)
            {
                if (rTreePanel1.VisibleLevels.Contains(e.Index + 1) == false)
                { 
                    rTreePanel1.VisibleLevels.Add(e.Index+1);
                    rTreePanel1.RefreshChart();
                }
            }
            else if (e.NewValue == CheckState.Unchecked)
            {
                if (rTreePanel1.VisibleLevels.Contains(e.Index + 1))
                { 
                    rTreePanel1.VisibleLevels.Remove(e.Index + 1);
                    rTreePanel1.RefreshChart();
                }
            }
        }

        private void checkedListBox2_ItemCheck(object sender, ItemCheckEventArgs e)
        {
            if (e.NewValue == CheckState.Checked)
            {
                if (rTreePanel1.VisibleGroups.Contains(e.Index + 1) == false)
                {
                    rTreePanel1.VisibleGroups.Add(e.Index + 1);
                    rTreePanel1.RefreshChart();
                }
            }
            else if (e.NewValue == CheckState.Unchecked)
            {
                if (rTreePanel1.VisibleGroups.Contains(e.Index + 1))
                {
                    rTreePanel1.VisibleGroups.Remove(e.Index + 1);
                    rTreePanel1.RefreshChart();
                }
            }
        }
    }
}
