using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Point = RTree.Point;

namespace mbR_NGQ
{
    public class UniformGenerator : DataGenerator
    {
        public override Point[] GeneratePoints(int pointNum)
        {
            Point[] points = new Point[pointNum];

            Random rand = new Random();

            for (int i = 0; i < pointNum; i++)
            {
                float x = (float) (rand.NextDouble() * (Config.maxX - Config.minX) + Config.minX);
                float y = (float) (rand.NextDouble() * (Config.maxY - Config.minY) + Config.minY);
                byte category = (byte)rand.Next(Config.m);

                points[i] = new Point(x, y, category);
            }

            return points;
        }
    }
}
