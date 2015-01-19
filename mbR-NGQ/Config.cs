using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mbR_NGQ
{
    public static class Config
    {
        public static int m = 3; // number of category
        public static int n = 100000; // number of data points, MAX 1M

        public static int k = 10; // number of Nearest Group

        public static float minX = 0;
        public static float maxX = 100;
        public static float minY = 0;
        public static float maxY = 100;
    }
}
