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
        public override Point[] GeneratePoints(int pointNum, int m, Random rand)
        {
            Point[] points = new Point[pointNum * m];

            for (int i = 0; i < pointNum * m; i++)
            {
                float x = (float)(rand.NextDouble() * (Config.maxX - Config.minX) + Config.minX);
                float y = (float)(rand.NextDouble() * (Config.maxY - Config.minY) + Config.minY);
                byte category = (byte)rand.Next(m);

                points[i] = new Point(x, y, category);
                points[i].id = i;
            }

            return points;
        }
    }

    public class GaussianIslandGenerator : DataGenerator
    {
        public override Point[] GeneratePoints(int pointNum, int m, Random rand)
        {
            Point[] points = new Point[pointNum * m];

            Point[] mPoints = new Point[m];
            for (int mIdx = 0; mIdx < m; mIdx++)
            {
                float x = (float)(rand.NextDouble() * (Config.maxX - Config.minX) * 0.8 + Config.minX * 0.2 * (Config.maxX - Config.minX)); // 10% margin
                float y = (float)(rand.NextDouble() * (Config.maxY - Config.minY) * 0.8 + Config.minY * 0.2 * (Config.maxX - Config.minX)); // 10% margin

                mPoints[mIdx] = new Point(x, y, 0);
            }

            int idx = 0;

            while (idx < pointNum * m)
            {
                byte category = (byte)rand.Next(m);

                float x = (float)NextGaussian(rand, (double)mPoints[category].coordinates[0], (Config.maxX - Config.minX) * 0.05);
                float y = (float)NextGaussian(rand, (double)mPoints[category].coordinates[1], (Config.maxY - Config.minY) * 0.05);

                if (x < Config.minX) continue;
                else if (x > Config.maxX) continue;
                else if (y < Config.minY) continue;
                else if (y > Config.maxY) continue;

                points[idx] = new Point(x, y, category);
                points[idx].id = idx;

                idx++;
            }

            return points;
        }

        //Box-Muller transform
        //code from : http://stackoverflow.com/questions/218060/random-gaussian-variables
        public double NextGaussian(Random r, double mu = 0, double sigma = 1)
        {
            var u1 = r.NextDouble();
            var u2 = r.NextDouble();

            var rand_std_normal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                                Math.Sin(2.0 * Math.PI * u2);

            var rand_normal = mu + sigma * rand_std_normal;

            return rand_normal;
        }

    }

    public class GaussianGenerator : DataGenerator
    {
        public override Point[] GeneratePoints(int pointNum, int m, Random rand)
        {
            Point[] points = new Point[pointNum * m];

            Point[] mPoints = new Point[m];
            for (int mIdx = 0; mIdx < m; mIdx++)
            {
                float x = (float)(rand.NextDouble() * (Config.maxX - Config.minX) * 0.2  + Config.minX * 0.4 * (Config.maxX - Config.minX)); // 40% margin
                float y = (float)(rand.NextDouble() * (Config.maxY - Config.minY) * 0.2 + Config.minY * 0.4 * (Config.maxX - Config.minX)); // 40% margin

                mPoints[mIdx] = new Point(x, y, 0);
            }

            int idx = 0;

            while (idx < pointNum * m)
            {
                byte category = (byte) rand.Next(m);

                float x = (float) NextGaussian(rand, (double)mPoints[category].coordinates[0], (Config.maxX - Config.minX) * 0.3);
                float y = (float)NextGaussian(rand, (double)mPoints[category].coordinates[1], (Config.maxY - Config.minY) * 0.3);

                if (x < Config.minX) continue;
                else if (x > Config.maxX) continue;
                else if (y < Config.minY) continue;
                else if (y > Config.maxY) continue;

                points[idx] = new Point(x, y, category);
                points[idx].id = idx;

                idx++;
            }

            return points;
        }

        //Box-Muller transform
        //code from : http://stackoverflow.com/questions/218060/random-gaussian-variables
        public double NextGaussian(Random r, double mu = 0, double sigma = 1)
        {
            var u1 = r.NextDouble();
            var u2 = r.NextDouble();

            var rand_std_normal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                                Math.Sin(2.0 * Math.PI * u2);

            var rand_normal = mu + sigma * rand_std_normal;

            return rand_normal;
        }

    }
}
