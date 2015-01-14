using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Point = RTree.Point;

namespace mbR_NGQ
{
    public abstract class DataGenerator
    {
        public abstract Point[] GeneratePoints(int n);
    }
}
