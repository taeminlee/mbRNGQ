//   Point.java
//   Java Spatial Index Library
//   Copyright (C) 2002 Infomatiq Limited
//  
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//  
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//  
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

// Ported to C# By Dror Gluska, April 9th, 2009

using System;
using System.Security.Cryptography.X509Certificates;

namespace RTree
{

    /// <summary>
    /// Currently hardcoded to 3 dimensions, but could be extended.
    /// author  aled@sourceforge.net
    /// version 1.0b2p1
    /// </summary>
    public class Point
    {
        public override string ToString()
        {
            return string.Format("ID : {0} x : {1} y : {2}", this.id, this.coordinates[0], this.coordinates[1]);
        }

        /// <summary>
        /// Number of dimensions in a point. In theory this
        /// could be exended to three or more dimensions.
        /// </summary>
        private const int DIMENSIONS = 2; // 3 to 2 ภฬลยนฮ 20150106

        internal byte category;

        /// <summary>
        /// The (x, y) coordinates of the point.
        /// </summary>
        internal float[] coordinates;

        public int id;

        /// <summary>
        /// Constructor.
        /// </summary>
        /// <param name="x">The x coordinate of the point</param>
        /// <param name="y">The y coordinate of the point</param>
        public Point(float x, float y, byte category)
        {
            coordinates = new float[DIMENSIONS];
            coordinates[0] = x;
            coordinates[1] = y;
            this.category = category;
        }

        public double GetDist(Point p)
        {
            double dist = 0;
            for (int i = 0; i < DIMENSIONS; i++)
            {
                dist += Math.Pow(p.coordinates[i] - this.coordinates[i],2);
            }
            return Math.Sqrt(dist);
        }
    }
}