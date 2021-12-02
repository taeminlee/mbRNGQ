using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using mbR_NGQ;
using RTree;

namespace UnitTest
{
    [TestClass]
    public class RTreeTest
    {
        [TestMethod]
        public void PointDist()
        {
            RTree.Point point1 = new Point(0, 0, 0);
            RTree.Point point2 = new Point(10, 10, 0);
            Assert.AreEqual(point1.GetDist(point2), Math.Sqrt(200), 0.0001);
        }

        public void MBRDist()
        {
            RTree.Rectangle r1 = new Rectangle(0,0, 10,10);
            RTree.Rectangle r2 = new Rectangle(0, 0, 10, 10);

            

            //Assert.AreEqual();
        }
    }
}
