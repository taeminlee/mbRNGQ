using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using RTree;
using Point = RTree.Point;
using System.Windows.Forms.DataVisualization.Charting;
using Rectangle = RTree.Rectangle;

namespace mbR_NGQ
{
    public partial class RTreePanel : UserControl
    {
        private RTree<Point> rTree;

        private List<int> visibleLevels = new List<int>();
        private List<int> visibleGroups = new List<int>();
        private List<int> visibleGoldGroups = new List<int>();
        private Point queryPoint = null;
        private bool showAllNodes = false;
        private List<RTree<Point>.MGroup<Point>> nearGroups = null;
        private List<RTree<Point>.MGroup<Point>> goldGroups = null;

        public List<int> VisibleLevels
        {
            get { return visibleLevels; }
        }

        public List<int> VisibleGroups
        {
            get { return visibleGroups; }
        }

        public List<int> VisibleGoldGroups
        {
            get { return visibleGoldGroups; }
        }
        
        public RTree<Point> RTree
        {
            get { return rTree; }
            set { rTree = value; }
        }

        public List<RTree<Point>.MGroup<Point>> NearGroups
        {
            get { return nearGroups; }
            set { nearGroups = value; }
        }

        public List<RTree<Point>.MGroup<Point>> GoldGroups
        {
            get { return goldGroups; }
            set { goldGroups = value; }
        }

        public Point QueryPoint
        {
            get { return queryPoint; }
            set { queryPoint = value; }
        }

        public bool ShowAllNodes
        {
            get { return showAllNodes;  }
            set { showAllNodes = value; }
        }

        public RTreePanel()
        {
            InitializeComponent();
        }

        private void InitChart()
        {
            ChartArea chartArea = this.chart1.ChartAreas.First();
            chartArea.AxisX.Minimum = Config.minX;
            chartArea.AxisX.Maximum = Config.maxX;
            chartArea.AxisY.Minimum = Config.minY;
            chartArea.AxisY.Maximum = Config.maxY;
        }

        public void RefreshChart()
        {
            this.chart1.Invalidate();
            this.chart1.Update();
            this.chart1.Refresh();
        }

        public void InitRTree()
        {
            if (rTree == null) return;

            InitChart();

            foreach (Point p in rTree.Contains(new Rectangle(Config.minX, Config.minY, Config.maxX, Config.maxY)))
            {
                int id = chart1.Series[0].Points.AddXY(p.coordinates[0], p.coordinates[1]);
                DataPoint dp = chart1.Series[0].Points[id];
                //dp.MarkerSize = 15;
                switch (p.category)
                {
                    case 0 :
                        dp.MarkerStyle = MarkerStyle.Circle;
                        dp.Color = Color.Blue;
                        break;
                    case 1:
                        dp.MarkerStyle = MarkerStyle.Cross;
                        dp.Color = Color.Gold;
                        break;
                    case 2:
                        dp.MarkerStyle = MarkerStyle.Diamond;
                        dp.Color = Color.Green;
                        break;
                    case 3:
                        dp.MarkerStyle = MarkerStyle.Square;
                        dp.Color = Color.Tan;
                        break;
                    case 4:
                        dp.MarkerStyle = MarkerStyle.Triangle;
                        dp.Color = Color.Orchid;
                        break;
                    case 5:
                        dp.MarkerStyle = MarkerStyle.Star10;
                        dp.Color = Color.PaleTurquoise;
                        break;
                }
            }

            if (queryPoint != null)
            {
                chart1.Series[1].Points.AddXY(queryPoint.coordinates[0], queryPoint.coordinates[1]);
            }
        }

        private void chart1_PostPaint(object sender, System.Windows.Forms.DataVisualization.Charting.ChartPaintEventArgs e)
        {
            if (rTree == null) return;

            //System.Diagnostics.Debug.WriteLine("CHART1 POSTPAINT");

            Graphics graphics = e.ChartGraphics.Graphics;
            foreach (var item in rTree.NodeMap)
            {
                Node<Point> node = item.Value;

                if (!showAllNodes && !node.opened) continue;

                if(visibleLevels.Contains(node.level))
                {
                    Pen pen = Pens.Black;
                    
                    switch (visibleLevels.IndexOf(node.level))
                    {
                        case 0 :
                            pen = Pens.CadetBlue;
                            break;
                        case 1:
                            pen = Pens.CornflowerBlue;
                            break;
                        case 2:
                            pen = Pens.BlueViolet;
                            break;
                        case 3:
                            pen = Pens.Blue;
                            break;
                        case 4:
                            pen = Pens.ForestGreen;
                            break;
                        case 5:
                            pen = Pens.Indigo;
                            break;
                    }
                    if (node.visited)
                    {
                        pen = Pens.Red;
                    }
                    if (!node.opened)
                    {
                        pen = Pens.Gray;
                    }
                    DrawMBR(e.ChartGraphics.Graphics, node.mbr, pen);
                    DrawBitArray(e.ChartGraphics.Graphics, node.mbr, node.bitArray);
                    DrawNodeID(e.ChartGraphics.Graphics, node.mbr, node.nodeId);
                }
            }

            chart1.Series[2].Points.Clear();
            foreach (int groupIdx in visibleGroups)
            {
                RTree<Point>.MGroup<Point> group = nearGroups[groupIdx - 1];
                foreach (Point point in group.Elems)
                {
                    chart1.Series[2].Points.AddXY(point.coordinates[0], point.coordinates[1]);
                }
            }

            foreach (int groupIdx in visibleGoldGroups)
            {
                RTree<Point>.MGroup<Point> group = goldGroups[groupIdx - 1];
                foreach (Point point in group.Elems)
                {
                    chart1.Series[2].Points.AddXY(point.coordinates[0], point.coordinates[1]);
                }
            }
        }

        private void DrawMBR(Graphics graphics, Rectangle mbr, Pen pen)
        {
            float x = (float)chart1.ChartAreas["rtree"].AxisX.ValueToPixelPosition(mbr.min[0]);
            float y = (float)chart1.ChartAreas["rtree"].AxisY.ValueToPixelPosition(mbr.max[1]);
            float width = (float)chart1.ChartAreas["rtree"].AxisX.ValueToPixelPosition(mbr.max[0]) - x;
            float height = (float)chart1.ChartAreas["rtree"].AxisY.ValueToPixelPosition(mbr.min[1]) - y;

            //System.Diagnostics.Debug.WriteLine("x : " +x + " y : " + y + " width : " + width + " height : " + height);

            graphics.DrawRectangle(pen, x, y, width, height);
        }

        private void DrawBitArray(Graphics graphics, Rectangle mbr, BitArray bitArray)
        {
            float x = (float)chart1.ChartAreas["rtree"].AxisX.ValueToPixelPosition(mbr.min[0]);
            float y = (float)chart1.ChartAreas["rtree"].AxisY.ValueToPixelPosition(mbr.max[1]);

            string mbStr = "";

            foreach (bool bit in bitArray)
            {
                if (bit)
                    mbStr += "1";
                else
                    mbStr += "0";
            }

            graphics.DrawString(mbStr, DefaultFont, Brushes.Black, x, y);
        }

        private void DrawNodeID(Graphics graphics, Rectangle mbr, int nodeId)
        {
            float x = (float)chart1.ChartAreas["rtree"].AxisX.ValueToPixelPosition(mbr.max[0]);
            float y = (float)chart1.ChartAreas["rtree"].AxisY.ValueToPixelPosition(mbr.max[1]);

            string mbStr = nodeId.ToString();

            graphics.DrawString(mbStr, DefaultFont, Brushes.Red, x, y);
        }
    }
}
