using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Point = RTree.Point;

namespace mbR_NGQ
{
    public class DataSet
    {
        private Point[] points;
        private DataGenerator dataGenerator;

        public Point[] Points
        {
            get { return points; }
            set { points = value; }
        }

        public DataGenerator DataGenerator
        {
            get { return dataGenerator; }
            set { dataGenerator = value; }
        }


        public void DataGeneration(int pointNum)
        {
            System.Diagnostics.Debug.Assert(pointNum > 0);

            this.points = dataGenerator.GeneratePoints(pointNum);
        }
    }
}
