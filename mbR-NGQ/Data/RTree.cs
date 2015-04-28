//   RTree.java
//   Java Spatial Index Library
//   Copyright (C) 2002 Infomatiq Limited
//   Copyright (C) 2008 Aled Morris aled@sourceforge.net
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
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

//  Ported to C# By Dror Gluska, April 9th, 2009


using System.Collections.Generic;
using System;
using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Drawing.Drawing2D;
using System.Linq;
using System.Net;
using System.Runtime.CompilerServices;
using System.Security.Permissions;
using System.Text.RegularExpressions;
using System.Windows.Forms;
using mbR_NGQ;
using mbR_NGQ.Data;
using mbR_NGQ.UI;
using Priority_Queue;
using Point = RTree.Point;


namespace RTree
{

    /// <summary>
    /// This is a lightweight RTree implementation, specifically designed 
    /// for the following features (in order of importance): 
    ///
    /// Fast intersection queryPoint performance. To achieve this, the RTree 
    /// uses only main memory to store entries. Obviously this will only improve
    /// performance if there is enough physical memory to avoid paging.
    /// Low memory requirements.
    /// Fast add performance.
    ///
    ///
    /// The main reason for the high speed of this RTree implementation is the 
    /// avoidance of the creation of unnecessary objects, mainly achieved by using
    /// primitive collections from the trove4j library.
    /// author aled@sourceforge.net
    /// version 1.0b2p1
    /// Ported to C# By Dror Gluska, April 9th, 2009
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class RTree<T>
    {

        private const string version = "1.0b2p1";

        // parameters of the tree
        private const int DEFAULT_MAX_NODE_ENTRIES = 10;
        internal int maxNodeEntries;
        int minNodeEntries;

        // map of nodeId -&gt; Node&lt;T&gt; object
        // [x] TODO eliminate this map - it should not be needed. Nodes
        // can be found by traversing the tree.
        //private TIntObjectHashMap nodeMap = new TIntObjectHashMap();
        private Dictionary<int, Node<T>> nodeMap = new Dictionary<int, Node<T>>();

        // to visualize mbr easily... , 20150105 이태민
        public Dictionary<int, Node<T>> NodeMap
        {
            get { return nodeMap; }
        }

        // internal consistency checking - set to true if debugging tree corruption
        private const bool INTERNAL_CONSISTENCY_CHECKING = false;

        // used to mark the status of entries during a Node&lt;T&gt; split
        private const int ENTRY_STATUS_ASSIGNED = 0;
        private const int ENTRY_STATUS_UNASSIGNED = 1;
        private byte[] entryStatus = null;
        private byte[] initialEntryStatus = null;

        // stacks used to store nodeId and entry index of each Node&lt;T&gt;
        // from the root down to the leaf. Enables fast lookup
        // of nodes when a split is propagated up the tree.
        //private TIntStack parents = new TIntStack();
        private Stack<int> parents = new Stack<int>();
        //private TIntStack parentsEntry = new TIntStack();
        private Stack<int> parentsEntry = new Stack<int>();

        // initialisation
        private int treeHeight = 1; // leaves are always level 1
        private int rootNodeId = 0;
        private int msize = 0;

        public int TreeHeight
        {
            get { return treeHeight; }
        }
        public int RootNodeId
        {
            get { return rootNodeId; }
        }

        // Enables creation of new nodes
        //private int highestUsedNodeId = rootNodeId; 
        private int highestUsedNodeId = 0;

        // Deleted Node&lt;T&gt; objects are retained in the nodeMap, 
        // so that they can be reused. Store the IDs of nodes
        // which can be reused.
        //private TIntStack deletedNodeIds = new TIntStack();
        private Stack<int> deletedNodeIds = new Stack<int>();

        // List of nearest rectangles. Use a member variable to
        // avoid recreating the object each time nearest() is called.
        //private TIntArrayList nearestIds = new TIntArrayList();
        List<int> nearestIds = new List<int>();

        //Added dictionaries to support generic objects..
        //possibility to change the code to support objects without dictionaries.
        private Dictionary<int, T> IdsToItems = new Dictionary<int,T>();
        private Dictionary<T,int> ItemsToIds = new Dictionary<T,int>();
        //private volatile int idcounter = int.MinValue;
        private volatile int idcounter = 0; // to know what is assigned easily... 20150106 이태민

        public Dictionary<int, T> IdsToItems1
        {
            get { return IdsToItems; }
        }

        public Dictionary<T, int> ItemsToIds1
        {
            get { return ItemsToIds; }
        }

        //the recursion methods require a delegate to retrieve data
        private delegate void intproc(int x);

        /// <summary>
        /// Initialize implementation dependent properties of the RTree.
        /// </summary>
        public RTree()
        {
            init();
        }

        /// <summary>
        /// Initialize implementation dependent properties of the RTree.
        /// </summary>
        /// <param name="MaxNodeEntries">his specifies the maximum number of entries
        ///in a openNode. The default value is 10, which is used if the property is
        ///not specified, or is less than 2.</param>
        /// <param name="MinNodeEntries">This specifies the minimum number of entries
        ///in a openNode. The default value is half of the MaxNodeEntries value (rounded
        ///down), which is used if the property is not specified or is less than 1.
        ///</param>
        public RTree(int MaxNodeEntries, int MinNodeEntries)
        {
            minNodeEntries = MinNodeEntries;
            maxNodeEntries = MaxNodeEntries;
            init();
        }

        private void init()
        {
            //initialize logs

            // Obviously a Node&lt;T&gt; with less than 2 entries cannot be split.
            // The Node&lt;T&gt; splitting algorithm will work with only 2 entries
            // per openNode, but will be inefficient.
            if (maxNodeEntries < 2)
            {
                maxNodeEntries = DEFAULT_MAX_NODE_ENTRIES;
            }

            // The MinNodeEntries must be less than or equal to (int) (MaxNodeEntries / 2)
            if (minNodeEntries < 1 || minNodeEntries > maxNodeEntries / 2)
            {
                minNodeEntries = maxNodeEntries / 2;
            }

            entryStatus = new byte[maxNodeEntries];
            initialEntryStatus = new byte[maxNodeEntries];

            for (int i = 0; i < maxNodeEntries; i++)
            {
                initialEntryStatus[i] = ENTRY_STATUS_UNASSIGNED;
            }

            Node<T> root = new Node<T>(rootNodeId, 1, maxNodeEntries);
            nodeMap.Add(rootNodeId, root);
        }

        public void CalculateBitArray()
        {
            Node<T> root = NodeMap[RootNodeId];
            root.bitArray = calculateBitArray(root);
        }

        private BitArray calculateBitArray(Node<T> node)
        {
            if(node.isLeaf() == false)
            { 
                BitArray[] nodeBitArrays = new BitArray[node.entryCount];
                for (int i = 0; i < node.entryCount; i++)
                {
                    nodeBitArrays[i] = calculateBitArray(NodeMap[node.ids[i]]);
                }
                BitArray bitArray = new BitArray(Config.m);
                foreach (BitArray nodeBitArray in nodeBitArrays)
                {
                    bitArray = bitArray.Or(nodeBitArray);
                }
                node.bitArray = bitArray;
                return bitArray;
            }
            else
            {
                BitArray bitArray = new BitArray(Config.m);
                for (int i = 0; i < node.entryCount; i++)
                {
                    var obj = IdsToItems[node.ids[i]];
                    if (obj.GetType() == typeof (RTree.Point))
                    {
                        RTree.Point p = obj as RTree.Point;
                        bitArray[p.category] = true;
                    }
                }
                node.bitArray = bitArray;
                return bitArray;
            }
        }

        /// <summary>
        /// Adds an item to the spatial index
        /// </summary>
        /// <param name="r"></param>
        /// <param name="item"></param>
        public void Add(Rectangle r, T item)
        {
            //System.Diagnostics.Debug.WriteLine("ADD : " + item.ToString());

            idcounter++;
            int id = idcounter;

            IdsToItems.Add(id, item);
            ItemsToIds.Add(item, id);

            add(r, id);
        }

        /// <summary>
        /// Adds an item to the spatial index
        /// </summary>
        /// <param name="r"></param>
        /// <param name="item"></param>
        public void Add(Rectangle r, T item, Point q)
        {
            //System.Diagnostics.Debug.WriteLine("ADD : " + item.ToString());

            r.minDist = r.distance(q);

            idcounter++;
            int id = idcounter;

            IdsToItems.Add(id, item);
            ItemsToIds.Add(item, id);

            add(r, id);
        }

        private void add(Rectangle r, int id)
        {
            
            add(r.copy(), id, 1);

            msize++;
        }

        /// <summary>
        /// Adds a new entry at a specified level in the tree
        /// </summary>
        /// <param name="r"></param>
        /// <param name="id"></param>
        /// <param name="level"></param>
        private void add(Rectangle r, int id, int level)
        {
            // I1 [Find position for new record] Invoke ChooseLeaf to select a 
            // leaf Node&lt;T&gt; L in which to place r
            Node<T> n = chooseNode(r, level);
            Node<T> newLeaf = null;

            // I2 [Add record to leaf openNode] If L has room for another entry, 
            // install E. Otherwise invoke SplitNode to obtain L and LL containing
            // E and all the old entries of L
            if (n.entryCount < maxNodeEntries)
            {
                n.addEntryNoCopy(r, id);
            }
            else
            {
                newLeaf = splitNode(n, r, id);
            }

            // I3 [Propagate changes upwards] Invoke AdjustTree on L, also passing LL
            // if a split was performed
            Node<T> newNode = adjustTree(n, newLeaf);

            // I4 [Grow tree taller] If Node&lt;T&gt; split propagation caused the root to 
            // split, create a new root whose children are the two resulting nodes.
            if (newNode != null)
            {
                int oldRootNodeId = rootNodeId;
                Node<T> oldRoot = getNode(oldRootNodeId);

                rootNodeId = getNextNodeId();
                treeHeight++;
                Node<T> root = new Node<T>(rootNodeId, treeHeight, maxNodeEntries);
                root.addEntry(newNode.mbr, newNode.nodeId);
                root.addEntry(oldRoot.mbr, oldRoot.nodeId);
                if (nodeMap.ContainsKey(rootNodeId) == false)
                {
                    nodeMap.Add(rootNodeId, root);
                }
                else
                {
                    nodeMap[rootNodeId] = root;
                }
            }

            if (INTERNAL_CONSISTENCY_CHECKING)
            {
                checkConsistency(rootNodeId, treeHeight, null);
            }
        }

        /// <summary>
        /// Deletes an item from the spatial index
        /// </summary>
        /// <param name="r"></param>
        /// <param name="item"></param>
        /// <returns></returns>
        public bool Delete(Rectangle r, T item)
        {
            //System.Diagnostics.Debug.WriteLine("DELETE : " + item.ToString());
            int id = ItemsToIds[item];

            bool success = delete(r, id);
            if (success == true)
            {
                IdsToItems.Remove(id);
                ItemsToIds.Remove(item);
            }
            return success;
        }

        private bool delete(Rectangle r, int id)
        {
            // FindLeaf algorithm inlined here. Note the "official" algorithm 
            // searches all overlapping entries. This seems inefficient to me, 
            // as an entry is only worth searching if it contains (NOT overlaps)
            // the rectangle we are searching for.
            //
            // Also the algorithm has been changed so that it is not recursive.

            // FL1 [Search subtrees] If root is not a leaf, check each entry 
            // to determine if it contains r. For each entry found, invoke
            // findLeaf on the Node&lt;T&gt; pointed to by the entry, until r is found or
            // all entries have been checked.
            parents.Clear();
            parents.Push(rootNodeId);

            parentsEntry.Clear();
            parentsEntry.Push(-1);
            Node<T> n = null;
            int foundIndex = -1;  // index of entry to be deleted in leaf

            while (foundIndex == -1 && parents.Count > 0)
            {
                n = getNode(parents.Peek());
                int startIndex = parentsEntry.Peek() + 1;

                if (!n.isLeaf())
                {
                    bool contains = false;
                    for (int i = startIndex; i < n.entryCount; i++)
                    {
                        if (n.entries[i].contains(r))
                        {
                            parents.Push(n.ids[i]);
                            parentsEntry.Pop();
                            parentsEntry.Push(i); // this becomes the start index when the child has been searched
                            parentsEntry.Push(-1);
                            contains = true;
                            break; // ie go to next iteration of while()
                        }
                    }
                    if (contains)
                    {
                        continue;
                    }
                }
                else
                {
                    foundIndex = n.findEntry(r, id);
                }

                parents.Pop();
                parentsEntry.Pop();
            } // while not found

            if (foundIndex != -1)
            {
                n.deleteEntry(foundIndex, minNodeEntries);
                condenseTree(n);
                msize--;
            }

            // shrink the tree if possible (i.e. if root Node&lt;T%gt; has exactly one entry,and that 
            // entry is not a leaf openNode, delete the root (it's entry becomes the new root)
            Node<T> root = getNode(rootNodeId);
            while (root.entryCount == 1 && treeHeight > 1)
            {
                root.entryCount = 0;
                rootNodeId = root.ids[0];
                treeHeight--;
                root = getNode(rootNodeId);
            }

            return (foundIndex != -1);
        }

        /// <summary>
        /// Retrieve nearest items to a point in radius furthestDistance
        /// </summary>
        /// <param name="p">Point of origin</param>
        /// <param name="furthestDistance">maximum distance</param>
        /// <returns>List of items</returns>
        public List<T> Nearest(Point p, float furthestDistance)
        {
            List<T> retval = new List<T>();
            nearest(p, delegate(int id)
            {
                retval.Add(IdsToItems[id]);
            }, furthestDistance);
            return retval;
        }


        private void nearest(Point p, intproc v, float furthestDistance)
        {
            Node<T> rootNode = getNode(rootNodeId);

            nearest(p, rootNode, furthestDistance);

            foreach (int id in nearestIds)
                v(id);
            nearestIds.Clear();
        }

        public class Dist : PriorityQueueNode
        {
            private double value;

            public double Value
            {
                get { return value; }
                set { this.value = value; }
            }

            public Dist(double value)
            {
                this.value = value;
            }
        }

        public class MGroup<T> : PriorityQueueNode, IComparable
        {
            private object[] _elems;
            private double minDist;
            private double maxDist;
            internal Node<T> masterNode;
            public bool isMix = false;

            public override string ToString()
            {
                string str = "Leaf : " + IsLeaf().ToString();
                foreach (object obj in this._elems)
                {
                    if (obj.GetType() == typeof (Node<T>))
                    {
                        Node<T> node = obj as Node<T>;

                        str += " " + node.nodeId;
                    }
                }
                str += " MinDist : " + this.minDist;
                str += " isMix : " + this.isMix;
                return str;
            }

            public int CompareTo(object obj)
            {
                MGroup<T> group1 = this;
                MGroup<T> group2 = (MGroup<T>) obj;
                if (group1.minDist > group2.minDist)
                    return 1;
                else if (group1.minDist < group2.minDist)
                    return -1;
                else
                    return 0;
            }

            public object[] Elems
            {
                get { return _elems; }
                set { _elems = value; }
            }

            public double MinDist
            {
                get { return minDist; }
            }

            public double MaxDist
            {
                get { return maxDist; }
            }

            public int Count()
            {
                int cnt = 0;
                foreach (object obj in _elems)
                {
                    if (obj != null)
                        cnt++;
                }
                return cnt;
            }

            public int NullCount()
            {
                int cnt = 0;
                foreach (object obj in _elems)
                {
                    if (obj == null)
                        cnt++;
                }
                return cnt;
            }

            public int GetFirstNullIndex()
            {
                int idx = 0;
                foreach (object obj in _elems)
                {
                    if (obj != null)
                        idx++;
                    else
                    {
                        return idx;
                    }
                }
                return -1;
            }

            public bool IsDictator()
            {
                object lastObj = null;
                foreach (object obj in this._elems)
                {
                    if (obj == null) return false;
                    if (lastObj == null) lastObj = obj;
                    else
                    {
                        if (lastObj != obj) return false;
                    }
                }
                return true;
            }

            /// <summary>
            /// Get a minkovski sum using theta margin
            /// </summary>
            /// <param name="theta"></param>
            /// <param name="queryPoint"></param>
            /// <returns></returns>
            public Rectangle GetCombinationMBR(double theta, Point queryPoint)
            {
                
                Rectangle mbr = new Rectangle(Config.minX, Config.minY, Config.maxX, Config.maxY);
                mbr = UnionLocalThetaMBR(mbr, queryPoint, theta);

                foreach (Object obj in this._elems)
                {
                    if (obj == null)
                        continue;
                    else
                    {
                        if (obj.GetType() == typeof (Point))
                        {
                            mbr = UnionLocalThetaMBR(mbr, (Point)obj, theta);
                        }
                        else if (obj.GetType() == typeof (Node<T>))
                        {
                            mbr = UnionLocalThetaMBR(mbr, (Node<T>)obj, theta);
                        }
                    }
                }
                
                return mbr;
            }


            private Rectangle UnionLocalThetaMBR(Rectangle mbr, Point point, double theta)
            {
                mbr.min[0] = Math.Max(mbr.min[0], point.coordinates[0] - (float)theta);
                mbr.min[1] = Math.Max(mbr.min[1], point.coordinates[1] - (float)theta);
                mbr.max[0] = Math.Min(mbr.max[0], point.coordinates[0] + (float)theta);
                mbr.max[1] = Math.Min(mbr.max[1], point.coordinates[1] + (float)theta);

                return mbr;
            }

            private Rectangle UnionLocalThetaMBR(Rectangle mbr, Node<T> node, double theta)
            {
                mbr.min[0] = Math.Max(mbr.min[0], node.mbr.min[0] - (float)theta);
                mbr.min[1] = Math.Max(mbr.min[1], node.mbr.min[1] - (float)theta);
                mbr.max[0] = Math.Min(mbr.max[0], node.mbr.max[0] + (float)theta);
                mbr.max[1] = Math.Min(mbr.max[1], node.mbr.max[1] + (float)theta);

                return mbr;
            }

            public MGroup(object[] _elems)
            {
                System.Diagnostics.Debug.Assert(_elems.Length == Config.m);
                this._elems = _elems;
            }

            public MGroup(Node<T> node)
            {
                this.masterNode = node;
                this._elems = new object[node.bitArray.Length];
                for (int bIdx = 0; bIdx < node.bitArray.Length; bIdx++)
                {
                    if(node.bitArray[bIdx])
                        _elems[bIdx] = node;
                }
            }

            public MGroup(Node<T> node, int mIdx)
            {
                this.masterNode = node;
                this._elems = new object[node.bitArray.Length];
                _elems[mIdx] = node;
            }

            public bool IsLeaf()
            {
                foreach (object obj in _elems)
                {
                    if (obj == null) return false;
                    if (obj.GetType() != typeof(Point))
                    {
                        return false;
                    }
                }
                return true;
            }

            public bool IsInternal()
            {
                foreach (object obj in _elems)
                {
                    if (obj == null) return false;
                    if (obj.GetType() != typeof(Node<T>))
                    {
                        return false;
                    }
                }
                return true;
            }

            public bool IsAllSet()
            {
                bool noneSet = false;
                BitArray ba = new BitArray(Config.m);

                foreach (Object elem in _elems)
                {
                    if (elem == null)
                    {
                        noneSet = true;
                    }
                    else if (elem.GetType() == typeof(Node<T>))
                    {
                        Node<T> node = elem as Node<T>;
                        ba.Or(node.bitArray);
                    }
                    else if (elem.GetType() == typeof(Point))
                    {
                        Point point = elem as Point;
                        ba[point.category] = true;
                    }
                }

                if (noneSet == false)
                {
                    foreach (bool b in ba)
                    {
                        System.Diagnostics.Debug.Assert(b);
                    }
                }

                foreach(object obj in _elems)
                    if (obj == null) return false;
                return true;
            }

            public bool IsBlank()
            {
                foreach (object obj in _elems)
                {
                    if (obj != null) return false;
                }
                return true;
            }

            public bool IsAllSame()
            {
                object temp_obj = null;
                foreach (object obj in _elems)
                    if (temp_obj == null) temp_obj = obj;
                    else if (temp_obj != obj) return false;
                return true;
            }

            public double GetGroupMaxDist(Point q)
            {
                if (IsLeaf() == false)
                {
                    Rectangle r = GetGroupMBR();
                    double innerDist = GetDiagonalDist(r);
                    //double interDist = GetMaxMinInterDist(this, q);
                    //double interDist = GetMaxMinInterDist(r, q);
                    double interDist = Double.MaxValue;
                    foreach (var elem in _elems)
                    {
                        double _interDist = GetMaxInterDist(elem, q);
                        if (interDist > _interDist)
                            interDist = _interDist;
                    }
                    //double interDist = GetMaxInterDist(r, q);
                    this.maxDist = innerDist + interDist;
                    return this.maxDist;
                }
                else // isLeaf() == true
                {
                    if (minDist != 0) this.maxDist = this.minDist;
                    else
                    {
                        this.maxDist = GetGroupMinDist(q);
                    }
                    return this.maxDist;
                }
            }

            public double GetMaxMinInterDist(MGroup<T> group, Point q)
            {
                double interDist = Double.MaxValue;

                foreach (object elem in group.Elems)
                {
                    if (elem.GetType() == typeof (Node<T>))
                    {
                        Node<T> node = elem as Node<T>;
                        double maxMinDist = GetMaxMinInterDist(node.mbr, q);
                        //double maxMinDist = GetMaxInterDist(openNode.mbr, q);
                        if (maxMinDist < interDist)
                            interDist = maxMinDist;
                    }
                    else if (elem.GetType() == typeof (Point))
                    {
                        Point point = elem as Point;
                        double maxMinDist = point.GetDist(q);
                        if (maxMinDist < interDist)
                            interDist = maxMinDist;
                    }
                }
                
                return interDist;
            }

            private double GetMaxInterDist(object elem, Point point)
            {
                if (elem.GetType() == typeof (Point))
                    return GetMaxInterDist((Point) elem, point);
                else if (elem.GetType() == typeof(Node<T>))
                    return GetMaxInterDist(((Node<T>)elem).mbr, point);
                System.Diagnostics.Debug.Fail("GetMaxInterDist Error");
                return -1;
            }

            // get maximal distance from query point to MBR, which is the maximum dist of interDist(q, M);
            public double GetMaxInterDist(Rectangle rectangle, Point point)
            {
                double[] center = new double[rectangle.min.Length];
                int dimension = center.Length;
                Point candidatePoint = new Point(0, 0, 0);
                for (int i = 0; i < dimension; i++)
                {
                    center[i] = (rectangle.min[i] + rectangle.max[i]) / 2;
                    if (point.coordinates[i] <= center[i])
                    {
                        candidatePoint.coordinates[i] = rectangle.max[i];
                    }
                    else
                    {
                        candidatePoint.coordinates[i] = rectangle.min[i];
                    }
                }
                double interDist = point.GetDist(candidatePoint);
                return interDist;
            }

            public double GetMaxInterDist(Point candidatePoint, Point point)
            {
                double interDist = point.GetDist(candidatePoint);
                return interDist;
            }

            // Find endpoints of Maximum digonal line, and then find minimum interDist of endpoints
            public double GetMaxMinInterDist(Rectangle rectangle, Point point)
            {
                double interDist = Double.MaxValue;
                double[] center = new double[rectangle.min.Length];
                int dimension = center.Length;
                BitArray pointSide = new BitArray(dimension);
                pointSide.SetAll(false);
                for (int i = 0; i < dimension; i++)
                {
                    center[i] = (rectangle.min[i] + rectangle.max[i]) / 2;
                    if (point.coordinates[i] > center[i]) pointSide[i] = true;
                    else pointSide[i] = false;
                }
                for (int i = 0; i < dimension; i++)
                {
                    BitArray candidatePointSide = new BitArray(dimension);
                    candidatePointSide.SetAll(false);
                    candidatePointSide[i] = true;
                    candidatePointSide = candidatePointSide.Xor(pointSide);
                    Point candidatePoint = new Point(0, 0, 0);
                    for (int j = 0; j < dimension; j++)
                    {
                        if (candidatePointSide[j] == true)
                        {
                            candidatePoint.coordinates[j] = rectangle.max[j];
                        }
                        else
                        {
                            candidatePoint.coordinates[j] = rectangle.min[j];
                        }
                    }
                    if (interDist > getMinDist(point, candidatePoint))
                    {
                        interDist = getMinDist(point, candidatePoint);
                    }
                }
                return interDist;
            }

            public double GetDiagonalDist(Rectangle rectangle)
            {
                double dist = 0;
                for (int i = 0; i < rectangle.min.Length; i++)
                {
                    dist += Math.Pow(rectangle.max[i] - rectangle.min[i], 2);
                }
                return Math.Sqrt(dist);
            }

            public Rectangle GetGroupMBR()
            {
                Rectangle r = null;

                if (IsAllSame()) return ((Node<T>)_elems[0]).mbr;

                foreach (var obj in _elems)
                {
                    if (obj.GetType() == typeof(Node<T>))
                    {
                        Node<T> node = (Node<T>)obj;
                        if (r == null)
                        {
                            r = node.mbr.copy();
                        }
                        else
                        {
                            r.add(node.mbr);
                        }
                    }
                    else if (obj.GetType() == typeof(Point))
                    {
                        Point point = (Point)obj;
                        if (r == null)
                        {
                            r = new Rectangle(point.coordinates, point.coordinates);
                        }
                        else
                        {
                            r.add(new Rectangle(point.coordinates, point.coordinates));
                        }
                    }
                    else
                    {
                        throw new NotImplementedException();
                    }
                }
                return r;
            }

            public double GetGroupMinDist(Point q)
            {
                double innerDist = 0;
                double interDist = Config.maxX - Config.minX + Config.maxY - Config.minY;
                for (int i = 0; i < _elems.Length; i++)
                {
                    if (_elems[i] == null) continue;

                    for (int j = i + 1; j < _elems.Length; j++)
                    {
                        if (_elems[i] == _elems[j]) continue;
                        if (_elems[j] == null) continue;

                        double ooDist = GetMinDist(_elems[i], _elems[j]);
                        if (innerDist < ooDist)
                        {
                            innerDist = ooDist;
                        }
                    }

                    double qoDist = GetMinDist(q, _elems[i]);
                    if (interDist > qoDist)
                    {
                        interDist = qoDist;
                    }
                }
                this.minDist = innerDist + interDist;
                return this.minDist;
            }

            private double GetMinDist(object obj1, object obj2)
            {
                if (obj1.GetType() == typeof(Node<T>) && obj2.GetType() == typeof(Node<T>))
                {
                    return getMinDist((Node<T>)obj1, (Node<T>)obj2);
                }
                else if (obj1.GetType() == typeof(Point) && obj2.GetType() == typeof(Point))
                {
                    return getMinDist((Point)obj1, (Point)obj2);
                }
                else if (obj1.GetType() == typeof(Node<T>) && obj2.GetType() == typeof(Point))
                {
                    return getMinDist((Node<T>)obj1, (Point)obj2);
                }
                else if (obj1.GetType() == typeof(Point) && obj2.GetType() == typeof(Node<T>))
                {
                    return getMinDist((Node<T>)obj2, (Point)obj1);
                }
                throw new NotImplementedException();
            }

            private double getMinDist(Node<T> obj1, Node<T> obj2)
            {
                return obj1.mbr.distance(obj2.mbr);
            }

            private double getMinDist(Node<T> obj1, Point obj2)
            {
                return obj1.mbr.distance(obj2);
            }

            private double getMinDist(Point obj1, Point obj2)
            {
                return obj1.GetDist(obj2);
            }

            public RTree<T>.MGroup<T> Copy()
            {
                RTree<T>.MGroup<T> newMGroup = new RTree<T>.MGroup<T>(new object[Config.m]);
                for (int i = 0; i < this._elems.Length; i++)
                {
                    newMGroup._elems[i] = this._elems[i];
                }
                newMGroup.isMix = this.isMix;
                newMGroup.masterNode = this.masterNode;
                return newMGroup;
            }
        }

        private double globalMaxDist;
        private double theta;
        private int KCounter = 0;
        //private List<MGroup<T>> candidateGroups;
        private HeapPriorityQueue<MGroup<T>> candidateGroups;
        //private Dictionary<object, List<MGroup<T>>> invertedCandidateGroups;
        //private Dictionary<object, List<MGroup<T>>> invertedCandidateGroups;
        private Dictionary<Node<T>, List<MGroup<T>>> invertedCandidateGroups;
        private MGroup<T>[] resultGroups;
        private RTree<object>[] invertedElemRTrees;
        private RTree<Point>[] invertedObjRTrees;
        private double[] thetaPool;
        private MGroup<T>[] thetaGroups;
        private int leafCount;
        private BitArray groupBitArray;
        private int accessCount;

        public List<MGroup<T>> NearestGroup(Point q)
        {
            //foreach()
            foreach (var nm in this.nodeMap)
            {
                nm.Value.visited = false;
                nm.Value.opened = false;
            }

            //System.Diagnostics.Debug.WriteLine("NEAREST GROUP QUERY START : " + q.ToString());

            List<MGroup<T>> retval = new List<MGroup<T>>();

            Node<T> rootNode = getNode(rootNodeId);
            rootNode.opened = true;

            globalMaxDist = new Point(Config.minX, Config.minY, 0).GetDist(new Point(Config.maxX, Config.maxY, 0));
            theta = globalMaxDist;
            thetaPool = new double[Config.k+1];
            thetaGroups = new MGroup<T>[Config.k+1];
            //candidateGroups = new List<MGroup<T>>();
            candidateGroups = new HeapPriorityQueue<MGroup<T>>(10000000);
            //invertedCandidateGroups = new Dictionary<object, List<MGroup<T>>>();
            invertedCandidateGroups = new Dictionary<Node<T>, List<MGroup<T>>>();

            resultGroups = new MGroup<T>[Config.k];
            // 차원별로 역인덱스 구축
            // Data RTree
            // Node RTree
            // Object RTree
            // 왜이렇게 겹치는게 많냐? - -;;;
            invertedElemRTrees = new RTree<object>[Config.m];
            invertedObjRTrees = new RTree<Point>[Config.m];
            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                invertedElemRTrees[mIdx] = new RTree<object>();
                invertedObjRTrees[mIdx] = new RTree<Point>();
            }
            KCounter = 0;
            for (int i = 0; i < thetaPool.Length; i++)
            {
                thetaPool[i] = globalMaxDist;
            }

            leafCount = 0;
            accessCount = 0;
            groupBitArray = new BitArray(Config.m);
            groupBitArray.SetAll(true);

            // 데이터 구조 초기화

            // 루트 노드 데이터 구축

            bool allSet = true;
            MGroup<T> rootGroup = new MGroup<T>(new object[Config.m]);
            rootGroup.masterNode = rootNode;
            Rectangle rootRec = new Rectangle(rootNode as Node<Point>);
            for (int idx = 0; idx < Config.m; idx++)
            {
                rootGroup.Elems[idx] = rootNode;
                invertedElemRTrees[idx].Add(rootRec, rootNode);
            }
            candidateGroups.Enqueue(rootGroup, rootGroup.GetGroupMinDist(q));
            invertedCandidateGroups[rootNode] = new List<MGroup<T>>();
            invertedCandidateGroups[rootNode].Add(rootGroup);

            // 준비 작업 완료

            // 실제 질의 처리 시작
            NearestGroup(q, new Node<T>[]{rootNode});

            // 결과 처리
            var v = resultGroups.ToList();
            v.Sort();
            resultGroups = v.ToArray();

            foreach (MGroup<T> ng in resultGroups)
            {
                retval.Add(ng);
            }

            return retval;
        }

        List<Node<T>> delNodes = new List<Node<T>>(); 
        List<MGroup<T>> delGroups = new List<MGroup<T>>();

        List<MGroup<T>> delCandidates = new List<MGroup<T>>();
        List<MGroup<T>> inCandidates = new List<MGroup<T>>();

        private void NearestGroup(Point q, Node<T>[] openNodes)
        {
            // 질의 점 q, 탐색 노드 openNode
            // 사용 데이터 구조 : candidate 저장 구조체
            //                  역 인덱스 (elem, openNode, obj)
            //                  theta pool, threshold...
            // 기본 idea : node를 열어서 새로운 조합 생성, 새로운 조합 마다 min dist 계산, min dist를 이용해서 threshold 정의, threshold 값을 이용해서 prunning, best openNode 선택, 재 탐색
            //System.Diagnostics.Debug.WriteLine("NearestGroup {0} inside openNode : {1}", accessCount, openNode.ToString());
            System.Diagnostics.Debug.WriteLine("theta {0} number of thetapool : {1}", theta, thetaPool.Count());
            //int cIdx = candidateGroups.Count > Config.k - (KCounter + 1) ? Config.k - (KCounter + 1) : candidateGroups.Count - 1;
            //candidateGroups[]
            System.Diagnostics.Debug.WriteLine("number of candidate groups : {0}", candidateGroups.Count);
            for(int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                Debug.WriteLine("Number of InvertedObj[{1}] Objects : {0}", invertedObjRTrees[mIdx].Count, mIdx);
                Debug.WriteLine("Number of InvertedElem[{1}] Objects : {0}", invertedElemRTrees[mIdx].Count, mIdx);
            }

            accessCount++;
            foreach(Node<T> openNode in openNodes)
            {
                Debug.Assert(!openNode.visited);
                openNode.visited = true;
                delNodes.Add(openNode);
            }

            //System.Diagnostics.Debug.WriteLine("Inverted Elem RTrees");
            //PrintRTree(invertedElemRTrees);
            //System.Diagnostics.Debug.WriteLine("Inverted Node RTrees");
            //PrintRTree(invertedNodeRTrees);

            //PrintCandidateSet();

            List<MGroup<T>> _beforecandidateGroups = candidateGroups.ToList();;

            Dictionary<Node<T>, List<MGroup<T>>> _beforeInvertedCandidateGroups = new Dictionary<Node<T>, List<MGroup<T>>>();
            foreach (var ttt in invertedCandidateGroups)
            {
                _beforeInvertedCandidateGroups.Add(ttt.Key, ttt.Value.ToList());
            }

            foreach (Node<T> openNode in openNodes)
            {
                Rectangle rec = new Rectangle(openNode as Node<Point>);

                // 1. openNode 들을 역 인덱스에서 제거
                for (int mIdx = 0; mIdx < Config.m; mIdx++)
                {
                    if (openNode.bitArray[mIdx] == true)
                    {
                        Debug.Assert(invertedElemRTrees[mIdx].ItemsToIds.ContainsKey(openNode));
                        invertedElemRTrees[mIdx].Delete(rec, openNode);
                    }
                }
            }
            
            List<object> newElems = new List<object>();
            List<object>[] invertedNewElems = new List<object>[Config.m];
            List<Node<T>> newNodes = new List<Node<T>>();
            List<Node<T>>[] invertedNewNodes = new List<Node<T>>[Config.m];
            List<Node<T>> dirtyNodes = new List<Node<T>>();

            foreach (Node<T> openNode in openNodes)
            {
                List<Point>[] invertedNewObjs = new List<Point>[Config.m];
                newElems = new List<object>();

                //4. Initialize Inverted List
                for (int mIdx = 0; mIdx < Config.m; mIdx++)
                {
                    invertedNewElems[mIdx] = new List<object>();
                    invertedNewObjs[mIdx] = new List<Point>();
                    invertedNewNodes[mIdx] = new List<Node<T>>();
                }

                System.Diagnostics.Debug.WriteLine("NearestGroup {0} inside openNode : {1}", accessCount,
                    openNode.ToString());

                List<KeyValuePair<Node<T>, MGroup<T>>> invDelGroups = new List<KeyValuePair<Node<T>, MGroup<T>>>();

                dirtyNodes = GetDirtyNodes(openNodes, openNode);

                // 5. openNode의 자식 object들을 메모리에 load
                GetOpenElements(q, openNode, invertedNewElems, newElems, invertedNewObjs, newNodes,
                    invertedNewNodes, openNodes);

                // 7. openNode의 자식 object와 기존 object와의 combination group을 생성

                // Make Combination Group (Computation Cost)

                if (openNode.isLeaf() == true) // Leaf Node
                {
                    MakeObjectGroups(q, invertedNewObjs, invertedObjRTrees, openNode);
                }

                // 6. openNode의 자식 object들을 역인덱스로 구축
                InsertOpenElements(q, newElems);
            }

                
            MakeBestGroups(q, newNodes, invertedNewNodes, invertedElemRTrees);

            // DirtyNode에 대해 BestGroup 생성 및 Theta Filtering 과정 진행. DirtyNode가 생성될 경우 해당 과정 반복
            while (dirtyNodes.Count > 0)
            {
                List<Node<T>>[] invertedDirtyNode = new List<Node<T>>[Config.m];

                for(int mIdx = 0; mIdx < Config.m; mIdx++)
                {
                   invertedDirtyNode[mIdx] = new List<Node<T>>();
                    foreach (Node<T> node in dirtyNodes)
                    {
                        if(node.bitArray[mIdx] == true)
                            invertedDirtyNode[mIdx].Add(node);
                    }
                }

                MakeBestGroups(q, dirtyNodes, invertedDirtyNode, invertedElemRTrees);

                dirtyNodes.Clear();

                //dirtyNodes = RemoveGroupByTheta(q);
            }

            foreach (Node<T> openNode in openNodes)
            {
                //RemoveElemInInvertedCandidateGroups(openNode);
                invertedCandidateGroups.Remove(openNode);
            }

            int currentKCounter = KCounter;

            List<MGroup<T>> _candidateGroups = candidateGroups.ToList();

            Dictionary<Node<T>, List<MGroup<T>>> _invertedCandidateGroups = new Dictionary<Node<T>, List<MGroup<T>>>();
            foreach (var ttt in invertedCandidateGroups)
            {
                _invertedCandidateGroups.Add(ttt.Key, ttt.Value.ToList());
            }


            List<Node<T>> bestNodes = new List<Node<T>>(Config.k);
            //double bestDist = globalMaxDist;

            //List<MGroup<T>> resultAddedGroups = new List<MGroup<T>>();

            int candidateGroupsCnt = candidateGroups.Count;

            for (int idx = 0; idx < candidateGroupsCnt; idx++)
            {
                //MGroup<T> group = candidateGroups[idx];
                MGroup<T> group = candidateGroups.Dequeue();
                if (group.IsLeaf())
                {
                    resultGroups[KCounter++] = group;

                    if (KCounter == Config.k)
                        return;

                    DeQueueThetaPool();
                    //resultAddedGroups.Add(group);
                }
                else
                {
                    delGroups.Add(group);
                    /*bestNodes.Add(group.masterNode);
                    break;*/
                    //if (bestNodes.Count > 0) continue;
                    foreach (object obj in group.Elems)
                    {
                        if (obj.GetType() == typeof(Node<T>))
                        {
                            Node<T> innerNode = (Node<T>)obj;
                            if(bestNodes.Contains(innerNode) == false)
                                bestNodes.Add(@innerNode);
                        }
                    }
                    break;
                }
            }

            /*foreach (MGroup<T> group in resultAddedGroups)
            {
                candidateGroups.Remove(group);
            }*/

            //List<MGroup<T>> _candidateGroups = candidateGroups.ToList();
            double _theta = theta;
            List<double> _thetaPool = thetaPool.ToList();
            int _kCounter = KCounter;

            System.Diagnostics.Debug.Assert(bestNodes.Count != 0);
            //System.Diagnostics.Debug.Assert(invertedCandidateGroups.ContainsKey(bestNode) != false);

            NearestGroup(q, bestNodes.ToArray());
        }

        /// <summary>
        /// dirty Nodes 찾기
        /// dirty Node는 openNode에 영향을 받는 openGroup에서 masterNode 중 openNodes에 포함되지 않는 것임
        /// 의미적으로 openNode에 의해 BestGroup이 변경되어 MinDist가 감소하게 되는 node를 의미함.
        /// </summary>
        /// <param name="openNodes"></param>
        /// <param name="openNode"></param>
        /// <returns></returns>
        private List<Node<T>> GetDirtyNodes(Node<T>[] openNodes, Node<T> openNode)
        {
            List<KeyValuePair<Node<T>, MGroup<T>>> invDelGroups = new List<KeyValuePair<Node<T>, MGroup<T>>>();
            List<Node<T>> dirtyNodes = new List<Node<T>>();
            List<int> checkNode = new List<int>();

            // 2. 
            foreach (MGroup<T> openGroup in invertedCandidateGroups[openNode])
            {
                ResetThetaPool(openGroup);
                candidateGroups.Remove(openGroup);
                delCandidates.Add(openGroup);

                checkNode.Clear();

                foreach (object obj in openGroup.Elems)
                {
                    if (obj.GetType() == typeof (Node<T>))
                    {
                        Node<T> node = (Node<T>) obj;
                        if (checkNode.Contains(node.nodeId)) continue;
                        checkNode.Add(node.nodeId);
                        invDelGroups.Add(new KeyValuePair<Node<T>, MGroup<T>>(node, openGroup));
                    }
                }

                if (dirtyNodes.Contains(openGroup.masterNode)) continue;
                if (openNodes.Contains(openGroup.masterNode)) continue;

                dirtyNodes.Add(openGroup.masterNode);
            }

            foreach (var group in invDelGroups)
            {
                bool tf = invertedCandidateGroups[group.Key].Remove(group.Value);
                Debug.Assert(tf);
            }

            return dirtyNodes;
        }

        public bool CheckNearGroup(Point queryPoint, List<MGroup<T>> nearGroups)
        {
            Rectangle range = null;
            List<Point> objects = new List<Point>();
            Dictionary<int, List<Point>> invertedObjects = new Dictionary<int, List<Point>>();
            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                invertedObjects[mIdx] = new List<Point>();
            }
            foreach (MGroup<T> group in nearGroups)
            {
                if (range == null)
                    range = group.GetGroupMBR();
                else
                {
                    range.add(group.GetGroupMBR());
                }
            }
            foreach (T t in this.Contains(range))
            {
                System.Diagnostics.Debug.Assert(t.GetType() == typeof(Point));

                Point point = t as Point;
                objects.Add(point);
                invertedObjects[point.category].Add(point);
            }

            List<MGroup<T>> groups = new List<MGroup<T>>();
            List<MGroup<T>> partialGroups = new List<MGroup<T>>();

            MGroup<T> templateGroup = new MGroup<T>(new object[Config.m]);
            partialGroups.Add(templateGroup);

            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                if (templateGroup.Elems[mIdx] != null) continue;

                List<MGroup<T>> removePartialGroups = new List<MGroup<T>>();
                List<MGroup<T>> addPartialGroups = new List<MGroup<T>>();
                foreach (object obj in invertedObjects[mIdx])
                {
                    foreach (MGroup<T> group in partialGroups)
                    {
                        if (group.Elems[mIdx] != null) continue;
                        removePartialGroups.Add(group);

                        MGroup<T> newGroup = group.Copy();
                        if (obj.GetType() == typeof(Node<T>))
                        {
                            Node<T> node = obj as Node<T>;
                            for (int bIdx = mIdx; bIdx < Config.m; bIdx++)
                            {
                                if (node.bitArray[bIdx] == true)
                                    newGroup.Elems[bIdx] = node;
                            }
                        }
                        else if (obj.GetType() == typeof(Point))
                        {
                            Point point = obj as Point;
                            newGroup.Elems[mIdx] = point;
                        }
                        addPartialGroups.Add(newGroup);
                    }
                }
                foreach (MGroup<T> group in removePartialGroups)
                {
                    partialGroups.Remove(group);
                }
                foreach (MGroup<T> group in addPartialGroups)
                {
                    if (group.IsAllSet())
                    {
                        group.GetGroupMinDist(queryPoint);
                        groups.Add(group);
                    }
                    else
                        partialGroups.Add(group);
                }
            }

            groups.Sort();

            for (int kIdx = 0; kIdx < Config.k; kIdx++)
            {
                if (resultGroups[kIdx].MinDist == groups[kIdx].MinDist) continue;
                else
                {
                    goldGroups = groups;
                    while (true)
                    {
                        this.NearestGroup(queryPoint);
                    }
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// open node가 포함된 group을 제거
        /// </summary>
        /// <param name="openNode">openNode</param>
        private void RemoveNodeInInvertedCandidateGroups(Node<T> openNode, Node<T>[] openNodes)
        {
            //System.Diagnostics.Debug.WriteLine("REMOVE ELEM " + elem);

            Debug.Assert(invertedCandidateGroups.ContainsKey(openNode));

            Dictionary<Node<T>, List<MGroup<T>>> removeGroups = new Dictionary<Node<T>, List<MGroup<T>>>();

            foreach (MGroup<T> group in invertedCandidateGroups[openNode])
            {
                List<object> dupList = new List<object>();
                foreach (object gElem in group.Elems)
                {
                    if (gElem.GetType() != typeof (Node<T>))
                        continue;
                    Node<T> node = (Node<T>) gElem;

                    if (openNodes.Contains(node)) continue;

                    if (dupList.Contains(node) == false)
                        dupList.Add(node);
                    if (removeGroups.ContainsKey(node) == false)
                        removeGroups[node] = new List<MGroup<T>>();
                    removeGroups[node].Add(group);
                }
            }
            foreach (var kv in removeGroups)
            {
                foreach (MGroup<T> group in kv.Value)
                {
                    //System.Diagnostics.Debug.WriteLine("REMOVE GROUP : " + kv.Key + " , " + group);
                    if (invertedCandidateGroups.ContainsKey(kv.Key))
                        if (invertedCandidateGroups[kv.Key].Contains(group))
                            invertedCandidateGroups[kv.Key].Remove(group);
                }
                if (invertedCandidateGroups.ContainsKey(kv.Key))
                    if (invertedCandidateGroups[kv.Key].Count == 0)
                        invertedCandidateGroups.Remove(kv.Key);
            }
        }

        private void MakeBestGroups(Point queryPoint, List<Node<T>> newNodes, List<Node<T>>[] invertedNewNodes, RTree<object>[] invertedElemRTrees)
        {
            foreach (Node<T> node in newNodes)
            {
                if (node.isDictator())
                {
                    MGroup<T> group = new MGroup<T>(node);
                    AddGroupToCandidate(queryPoint, group);
                }
                else
                {
                    MakeBestGroups(queryPoint, node, invertedNewNodes, invertedElemRTrees);
                }
            }
        }

        private void MakeBestGroups(Point queryPoint, Node<T> node, List<Node<T>>[] invertedNewNodes, RTree<object>[] invertedElemRTrees)
        {
            Debug.WriteLine("Make Best Groups {0}", node);
            Stopwatch sw = new Stopwatch();
            sw.Start();

            List<MGroup<T>> partialGroups = new List<MGroup<T>>();

            MGroup<T> bestGroup = null;

            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                if(node.bitArray[mIdx])
                { 
                    MGroup<T> nodeGroup = new MGroup<T>(node, mIdx);
                    partialGroups.Add(nodeGroup);
                }
            }

            while (partialGroups.Count != 0)
            {
                List<MGroup<T>> removePartialGroups = new List<MGroup<T>>();
                List<MGroup<T>> addPartialGroups = new List<MGroup<T>>();

                foreach (MGroup<T> partialGroup in partialGroups)
                {
                    int mIdx = partialGroup.GetFirstNullIndex();

                    Rectangle searchArea = partialGroup.GetCombinationMBR(theta, queryPoint);

                    System.Diagnostics.Debug.Assert(mIdx != -1);

                    foreach (object element in invertedElemRTrees[mIdx].Contains(searchArea))
                    {
                        System.Diagnostics.Debug.Assert(element != null);
                        if (element == node) continue;

                        MGroup<T> newGroup = partialGroup.Copy();
                        AddElementToGroup(element, newGroup, mIdx);

                        if (newGroup.IsAllSet())
                        {
                            newGroup.GetGroupMinDist(queryPoint);
                            if (bestGroup == null || bestGroup.MinDist > newGroup.MinDist)
                            {
                                bestGroup = newGroup;
                            }
                        }
                        else
                        {
                            if (partialGroup.GetGroupMinDist(queryPoint) > theta) continue;
                            addPartialGroups.Add(newGroup);
                        }
                    }
                    removePartialGroups.Add(partialGroup);
                }

                foreach (MGroup<T> group in removePartialGroups)
                {
                    partialGroups.Remove(group);
                }
                foreach (MGroup<T> group in addPartialGroups)
                {
                    partialGroups.Add(group);
                }
            }

            if (bestGroup == null) return;

            Debug.Assert(bestGroup != null);
            
            foreach (object elem in bestGroup.Elems)
            {
                if (elem.GetType() == typeof (Node<T>))
                {
                    Node<T> n = (Node<T>) elem;
                    Debug.Assert(!n.visited);
                }
            }
            AddGroupToCandidate(queryPoint, bestGroup);
            

            sw.Stop();
            Debug.WriteLine("Make Best Group Elapsed Time : {0}", sw.ElapsedMilliseconds);
        }

        private void MakeObjectGroups(Point queryPoint, List<Point>[] NewObjs, RTree<Point>[] objRTrees, Node<T> openNode)
        {
            BitArray combineBitArray = new BitArray(Config.m);
            combineBitArray.SetAll(true);

            for (int i = 0; i < Config.m; i++)
            {
                foreach (
                    var combineBitRule in
                        Util.Combinations(Enumerable.Range(0, Config.m).ToArray(), Config.m - i))
                {
                    bool needToRun = true;
                    combineBitArray.SetAll(false);
                    foreach (var cIdx in combineBitRule)
                    {
                        combineBitArray[cIdx] = true;
                        if (openNode.bitArray[cIdx] == false) needToRun = false;
                    }
                    if(needToRun)
                    { 
                        MakeObjectGroups(queryPoint, NewObjs, objRTrees, combineBitArray);
                    }
                }
            }
        }

        private void MakeObjectGroups(Point queryPoint, List<Point>[] NewObjs, RTree<Point>[] objRTrees, BitArray combineBitArray)
        {
            List<MGroup<T>> partialGroups = new List<MGroup<T>>();
            
            int groupCnt = 0;
            bool isMix = false;

            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                if (combineBitArray[mIdx] == true)
                {
                    foreach (Point point in NewObjs[mIdx])
                    {
                        MGroup<T> newGroup = new MGroup<T>(new object[Config.m]);
                        AddElementToGroup(point, newGroup, mIdx);
                        partialGroups.Add(newGroup);
                        
                    }
                    break;
                }
            }

            while (partialGroups.Count != 0)
            {
                List<MGroup<T>> removePartialGroups = new List<MGroup<T>>();
                List<MGroup<T>> addPartialGroups = new List<MGroup<T>>();

                foreach (MGroup<T> partialGroup in partialGroups)
                {
                    if (partialGroup.IsAllSet())
                    {
                        if (AddGroupToCandidate(queryPoint, partialGroup))
                        {
                            groupCnt++;
                        }
                        removePartialGroups.Add(partialGroup);
                        continue;
                    }

                    int mIdx = partialGroup.GetFirstNullIndex();

                    if (combineBitArray[mIdx] == false)
                    {
                        isMix = true;

                        Rectangle searchArea = partialGroup.GetCombinationMBR(theta, queryPoint);

                        System.Diagnostics.Debug.Assert(mIdx != -1);

                        foreach (object element in objRTrees[mIdx].Contains(searchArea))
                        {
                            System.Diagnostics.Debug.Assert(element != null);

                            MGroup<T> newGroup = partialGroup.Copy();
                            newGroup.isMix = true;
                            AddElementToGroup(element, newGroup, mIdx);
                            addPartialGroups.Add(newGroup);
                        }
                        removePartialGroups.Add(partialGroup);
                    }
                    else
                    {
                        foreach (object element in NewObjs[mIdx])
                        {
                            System.Diagnostics.Debug.Assert(element != null);

                            MGroup<T> newGroup = partialGroup.Copy();
                            AddElementToGroup(element, newGroup, mIdx);
                            addPartialGroups.Add(newGroup);
                        }
                        removePartialGroups.Add(partialGroup);
                    }
                }

                foreach (MGroup<T> group in removePartialGroups)
                {
                    partialGroups.Remove(group);
                }
                foreach (MGroup<T> group in addPartialGroups)
                {
                    partialGroups.Add(group);
                }
            }
            
            Debug.WriteLine("Combination Object Group, Mix : {0}, Cnt : {1}", isMix, groupCnt);
        }

        private void InsertOpenElements(Point q, List<object> newElem)
        {
            foreach (object elem in newElem)
            {
                if (elem.GetType() == typeof(Node<T>))
                {
                    Node<T> node = elem as Node<T>;
                    Rectangle elemRec = new Rectangle((Node<Point>)elem);
                    for (int mIdx = 0; mIdx < Config.m; mIdx++)
                    {
                        if(node.bitArray[mIdx] == true)
                        { 
                            invertedElemRTrees[mIdx].Add(elemRec, elem, q);
                        }
                    }
                }
                else if (elem.GetType() == typeof(Point))
                {
                    Rectangle elemRec = new Rectangle((Point)elem);
                    {
                        Point point = elem as Point;
                        invertedElemRTrees[((Point)elem).category].Add(elemRec, elem, q);
                        invertedObjRTrees[point.category].Add(elemRec, point, q);
                    }
                }
            }
        }

        private void ResetThetaPool(MGroup<T> group)
        {
            if (thetaGroups.Contains(group))
            {
                int thetaIdx = thetaGroups.ToList().IndexOf(group);
                thetaPool[thetaIdx] = theta;
            }
        }

        private void DeQueueThetaPool()
        {
            List<double> thetaPoolList = thetaPool.ToList();
            List<MGroup<T>> thetaGroupList = thetaGroups.ToList();
            int thetaIdx = thetaPool.ToList().IndexOf(thetaPool.Min());
            thetaPoolList.RemoveAt(thetaIdx);
            thetaGroupList.RemoveAt(thetaIdx);
            thetaPool = thetaPoolList.ToArray();
            thetaGroups = thetaGroupList.ToArray();
        }

        private void GetOpenElements(Point q, Node<T> node, List<object>[] invertedNewElems, List<object> newElems, List<Point>[] invertedNewObjs, List<Node<T>> newNodes, List<Node<T>>[] invertedNewNodes, Node<T>[] openNodes)
        {
            // Insert child nodes into inverted List
            if (node.isLeaf() == false)
            {
                Node<T>[] nodes = new Node<T>[node.entryCount];
                for (int i = 0; i < node.entryCount; i++)
                {
                    nodes[i] = NodeMap[node.ids[i]];
                    nodes[i].opened = true;
                }
                foreach (Node<T> childNode in nodes)
                {
                    if (childNode.mbr.distance(q) > theta) continue;

                    if (openNodes.Contains(childNode)) continue;

                    for (int mIdx = 0; mIdx < Config.m; mIdx++)
                    {
                        if (childNode.bitArray[mIdx] == true)
                        {
                            invertedNewElems[mIdx].Add(childNode);
                            invertedNewNodes[mIdx].Add(childNode);
                        }
                    }
                    newElems.Add(childNode);
                    newNodes.Add(childNode);
                }
            }
            // Add points into inverted List
            else if (node.isLeaf() == true)
            {
                T[] objs = new T[node.entryCount];
                for (int i = 0; i < node.entryCount; i++)
                {
                    objs[i] = IdsToItems[node.ids[i]];
                }
                foreach (T obj in objs)
                {
                    System.Diagnostics.Debug.Assert(obj.GetType() == typeof(Point));

                    if (obj.GetType() != typeof(Point)) continue;

                    Point point = obj as Point;

                    if (point.GetDist(q) > theta) continue;

                    invertedNewElems[point.category].Add(point);
                    invertedNewObjs[point.category].Add(point);
                    newElems.Add(point);
                }
            }
        }

        private bool AddGroupToCandidate(Point queryPoint, MGroup<T> @group)
        {
            //if (DuplicateCheck(@group) == false) return false;

            @group.GetGroupMinDist(queryPoint);
            if (@group.MinDist > theta)
                return false;
            @group.GetGroupMaxDist(queryPoint);
            UpdateTheta(@group);
            candidateGroups.Enqueue(@group, @group.MinDist);
            inCandidates.Add(@group);
            foreach(object elem in @group.Elems)
            {
                /*if (invertedCandidateGroups.ContainsKey(elem) == false)
                    invertedCandidateGroups[elem] = new List<MGroup<T>>();
                if (invertedCandidateGroups[elem].Contains(@group) == false)
                    invertedCandidateGroups[elem].Add(@group);*/
                if (elem.GetType() == typeof (Node<T>))
                {
                    Node<T> node = elem as Node<T>;
                    if (invertedCandidateGroups.ContainsKey(node) == false)
                        invertedCandidateGroups[node] = new List<MGroup<T>>();
                    if (invertedCandidateGroups[node].Contains(@group) == false)
                        invertedCandidateGroups[node].Add(@group);
                }
            }

            /*// All Set Test
            BitArray ba = new BitArray(Config.m);

            foreach (Object elem in @group.Elems)
            {
                if (elem.GetType() == typeof(Node<T>))
                {
                    Node<T> node = elem as Node<T>;
                    ba.Or(node.bitArray);
                }
                else if (elem.GetType() == typeof(Point))
                {
                    Point point = elem as Point;
                    ba[point.category] = true;
                }
            }

            foreach (bool b in ba)
            {
                System.Diagnostics.Debug.Assert(b);
            }*/

            return true;
        }

        private bool AddGroupToTheta(Point queryPoint, MGroup<T> @group)
        {
            @group.GetGroupMinDist(queryPoint);
            if (@group.MinDist > theta)
                return false;
            @group.GetGroupMaxDist(queryPoint);
            UpdateTheta(@group);
            return true;
        }

        /*private bool AddGroupToTheta(Point queryPoint, MGroup<T> @group)
        {
            @group.GetGroupMinDist(queryPoint);
            if (@group.MinDist > theta)
                return false;
            @group.GetGroupMaxDist(queryPoint);
            UpdateTheta(@group);
            foreach (object elem in @group.Elems)
            {
                if (invertedCandidateGroups.ContainsKey(elem) == false)
                    invertedCandidateGroups[elem] = new List<MGroup<T>>();
                invertedCandidateGroups[elem].Add(@group);
            }
            return true;
        }*/

        private bool DuplicateCheck(MGroup<T> partialGroup)
        {
            foreach (MGroup<T> group in candidateGroups)
            {
                if (IsDuplicated(partialGroup, @group) == true) return false;
            }
            foreach (MGroup<T> group in resultGroups)
            {
                if (group == null) continue;
                if (IsDuplicated(partialGroup, @group) == true) return false;
            }
            return true;
        }

        private bool IsDuplicated(MGroup<T> partialGroup, MGroup<T> @group)
        {
            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                if (partialGroup.Elems[mIdx].GetType() == typeof (Node<T>))
                {
                    if (@group.Elems[mIdx].GetType() != typeof (Node<T>))
                    {
                        return false;
                    }
                    Node<T> node1 = partialGroup.Elems[mIdx] as Node<T>;
                    Node<T> node2 = @group.Elems[mIdx] as Node<T>;

                    if (node1.nodeId != node2.nodeId)
                    {
                        return false;
                    }
                }
                if (partialGroup.Elems[mIdx].GetType() == typeof (Point))
                {
                    if (@group.Elems[mIdx].GetType() != typeof (Point))
                    {
                        return false;
                    }
                    Point obj1 = partialGroup.Elems[mIdx] as Point;
                    Point obj2 = @group.Elems[mIdx] as Point;

                    if (obj1.id != obj2.id)
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        private void AddElementToGroup(object element, MGroup<T> @group, int mIdx)
        {
            @group.Elems[mIdx] = element;
        }

        /// <summary>
        /// Threshold 값 이하의 MinDist를 갖는 Group과 Object를 제거함
        /// </summary>
        /// <param name="q"></param>
        /// <returns>재 계산이 필요한 Node 집합 (dirtyNodes)</returns>
        private List<Node<T>> RemoveGroupByTheta(Point q)
        {
            System.Diagnostics.Debug.WriteLine("RemoveGroupByTheta Candidate Num : {0}", candidateGroups.Count);

            List<Node<T>> dirtyNodes = new List<Node<T>>();

            Stopwatch sw = new Stopwatch();
            sw.Start();

            List<MGroup<T>> removeGroups = new List<MGroup<T>>();
            List<int> removeElements = new List<int>();
            foreach (var group in candidateGroups)
            {
                if (group.MinDist > theta)
                {
                    removeGroups.Add(group);
                    foreach (object elem in group.Elems)
                    {
                        if (elem.GetType() == typeof(Node<T>))
                        {
                            Node<T> node = elem as Node<T>;
                            if (node.mbr.distance(q) > theta)
                            {
                                if (removeElements.Contains(node.nodeId)) continue;
                                if (removeElements.Contains(node.nodeId) == false)
                                    removeElements.Add(node.nodeId);

                                System.Diagnostics.Debug.WriteLine("Remove Node : " + node.ToString());

                                Rectangle elemRec = new Rectangle((Node<Point>)elem);
                                for (int mIdx = 0; mIdx < Config.m; mIdx++)
                                {
                                    if(node.bitArray[mIdx] == true)
                                    { 
                                        invertedElemRTrees[mIdx].Delete(elemRec, elem);
                                    }
                                }
                                if (invertedCandidateGroups.ContainsKey(node) == false) continue;
                                foreach (MGroup<T> candidateGroup in invertedCandidateGroups[node])
                                {
                                    removeGroups.Add(candidateGroup);
                                    if(candidateGroup.masterNode != node) dirtyNodes.Add(node);
                                }
                                
                                //RemoveElemInInvertedCandidateGroups(node);
                            }
                        }
                        else if (elem.GetType() == typeof(Point))
                        {
                            Point point = elem as Point;
                            if(point.GetDist(q) > theta)
                            {
                                if (removeElements.Contains(point.id)) continue;
                                if (removeElements.Contains(point.id) == false)
                                    removeElements.Add(point.id);

                                Rectangle elemRec = new Rectangle((Point)elem);
                                invertedElemRTrees[point.category].Delete(elemRec, elem);
                                invertedObjRTrees[point.category].Delete(elemRec, point);
                            }
                        }
                    }
                }
            }

            sw.Stop();

            Debug.WriteLine("Remove Theta ElapsedTime : {0}", sw.ElapsedMilliseconds);

            //System.Diagnostics.Debug.Assert(sw.ElapsedMilliseconds < 2);

            foreach (var group in removeGroups)
            {
                candidateGroups.Remove(group);
            }

            return dirtyNodes;
        }

        private void UpdateTheta(MGroup<T> newGroup)
        {
            if (thetaPool.Max() > newGroup.MaxDist)
            {
                int idx = thetaPool.ToList().IndexOf(thetaPool.Max());
                thetaPool[idx] = newGroup.MaxDist;
                theta = thetaPool.Max();
                thetaGroups[idx] = newGroup;
            }
        }

        public List<MGroup<T>> goldGroups;


        /*public bool CheckNearGroup(Point queryPoint, List<MGroup<T>> nearGroups)
        {
            Rectangle range = null;
            List<Point> objects = new List<Point>();
            Dictionary<int, List<Point>> invertedObjects = new Dictionary<int, List<Point>>();
            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                invertedObjects[mIdx] = new List<Point>();
            }
            foreach (MGroup<T> group in nearGroups)
            {
                if (range == null)
                    range = group.GetGroupMBR();
                else
                {
                    range.add(group.GetGroupMBR());
                }
            }
            foreach (T t in this.Contains(range))
            {
                System.Diagnostics.Debug.Assert(t.GetType() == typeof(Point));

                Point point = t as Point;
                objects.Add(point);
                invertedObjects[point.category].Add(point);
            }

            List<MGroup<T>> groups = new List<MGroup<T>>();
            List<MGroup<T>> partialGroups = new List<MGroup<T>>();

            MGroup<T> templateGroup = new MGroup<T>(new object[Config.m]);
            partialGroups.Add(templateGroup);

            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                if (templateGroup.Elems[mIdx] != null) continue;

                List<MGroup<T>> removePartialGroups = new List<MGroup<T>>();
                List<MGroup<T>> addPartialGroups = new List<MGroup<T>>();
                foreach (object obj in invertedObjects[mIdx])
                {
                    foreach (MGroup<T> group in partialGroups)
                    {
                        if (group.Elems[mIdx] != null) continue;
                        removePartialGroups.Add(group);

                        MGroup<T> newGroup = group.Copy();
                        if (obj.GetType() == typeof(Node<T>))
                        {
                            Node<T> openNode = obj as Node<T>;
                            for (int bIdx = mIdx; bIdx < Config.m; bIdx++)
                            {
                                if (openNode.bitArray[bIdx] == true)
                                    newGroup.Elems[bIdx] = openNode;
                            }
                        }
                        else if (obj.GetType() == typeof(Point))
                        {
                            Point point = obj as Point;
                            newGroup.Elems[mIdx] = point;
                        }
                        addPartialGroups.Add(newGroup);
                    }
                }
                foreach (MGroup<T> group in removePartialGroups)
                {
                    partialGroups.Remove(group);
                }
                foreach (MGroup<T> group in addPartialGroups)
                {
                    if (group.IsAllSet())
                    {
                        group.GetGroupMinDist(queryPoint);
                        groups.Add(group);
                    }
                    else
                        partialGroups.Add(group);
                }
            }

            groups.Sort();

            for (int kIdx = 0; kIdx < Config.k; kIdx++)
            {
                if (resultGroups[kIdx].MinDist == groups[kIdx].MinDist) continue;
                else
                {
                    goldGroups = groups;
                    this.NearestGroup(queryPoint);
                    return false;
                }
            }

            return true;
        }*/


        /// <summary>
        /// Retrieve items which intersect with Rectangle r
        /// </summary>
        /// <param name="r"></param>
        /// <returns></returns>
        public List<T> Intersects(Rectangle r)
        {
            List<T> retval = new List<T>();
            intersects(r, delegate(int id)
            {
                retval.Add(IdsToItems[id]);
            });
            return retval;
        }


        private void intersects(Rectangle r, intproc v)
        {
            Node<T> rootNode = getNode(rootNodeId);
            intersects(r, v, rootNode);
        }

        /// <summary>
        /// find all rectangles in the tree that are contained by the passed rectangle
        /// written to be non-recursive (should model other searches on this?)</summary>
        /// <param name="r"></param>
        /// <returns></returns>
        public List<T> Contains(Rectangle r)
        {
            List<T> retval = new List<T>();
            contains(r, delegate(int id)
            {
                retval.Add(IdsToItems[id]);
            });

            return retval;
        }

        private void contains(Rectangle r, intproc v)
        {
            // find all rectangles in the tree that are contained by the passed rectangle
            // written to be non-recursive (should model other searches on this?)

            parents.Clear();
            parents.Push(rootNodeId);

            parentsEntry.Clear();
            parentsEntry.Push(-1);

            // TODO: possible shortcut here - could test for intersection with the 
            // MBR of the root openNode. If no intersection, return immediately.

            while (parents.Count > 0)
            {
                Node<T> n = getNode(parents.Peek());
                int startIndex = parentsEntry.Peek() + 1;

                if (!n.isLeaf())
                {
                    // go through every entry in the index Node<T> to check
                    // if it intersects the passed rectangle. If so, it 
                    // could contain entries that are contained.
                    bool intersects = false;
                    for (int i = startIndex; i < n.entryCount; i++)
                    {
                        if (r.intersects(n.entries[i]))
                        {
                            parents.Push(n.ids[i]);
                            parentsEntry.Pop();
                            parentsEntry.Push(i); // this becomes the start index when the child has been searched
                            parentsEntry.Push(-1);
                            intersects = true;
                            break; // ie go to next iteration of while()
                        }
                    }
                    if (intersects)
                    {
                        continue;
                    }
                }
                else
                {
                    // go through every entry in the leaf to check if 
                    // it is contained by the passed rectangle
                    for (int i = 0; i < n.entryCount; i++)
                    {
                        if (r.contains(n.entries[i]))
                        {
                            v(n.ids[i]);
                        }
                    }
                }
                parents.Pop();
                parentsEntry.Pop();
            }
        }

        /**
        * @see com.infomatiq.jsi.SpatialIndex#getBounds()
        */
        public Rectangle getBounds()
        {
            Rectangle bounds = null;

            Node<T> n = getNode(getRootNodeId());
            if (n != null && n.getMBR() != null)
            {
                bounds = n.getMBR().copy();
            }
            return bounds;
        }

        /**
         * @see com.infomatiq.jsi.SpatialIndex#getVersion()
         */
        public string getVersion()
        {
            return "RTree-" + version;
        }
        //-------------------------------------------------------------------------
        // end of SpatialIndex methods
        //-------------------------------------------------------------------------

        /**
         * Get the next available Node&lt;T&gt; ID. Reuse deleted Node&lt;T&gt; IDs if
         * possible
         */
        private int getNextNodeId()
        {
            int nextNodeId = 0;
            if (deletedNodeIds.Count > 0)
            {
                nextNodeId = deletedNodeIds.Pop();
            }
            else
            {
                nextNodeId = 1 + highestUsedNodeId++;
            }
            return nextNodeId;
        }

       
       
       

        /// <summary>
        /// Get a Node&lt;T&gt; object, given the ID of the openNode.
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        private Node<T> getNode(int index)
        {
            return (Node<T>)nodeMap[index];
        }

        /// <summary>
        /// Get the highest used Node&lt;T&gt; ID
        /// </summary>
        /// <returns></returns>
        private int getHighestUsedNodeId()
        {
            return highestUsedNodeId;
        }

        /// <summary>
        /// Get the root Node&lt;T&gt; ID
        /// </summary>
        /// <returns></returns>
        public int getRootNodeId()
        {
            return rootNodeId;
        }

        /// <summary>
        /// Split a openNode. Algorithm is taken pretty much verbatim from
        /// Guttman's original paper.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="newRect"></param>
        /// <param name="newId"></param>
        /// <returns>return new Node&lt;T&gt; object.</returns>
        private Node<T> splitNode(Node<T> n, Rectangle newRect, int newId)
        {
            // [Pick first entry for each group] Apply algorithm pickSeeds to 
            // choose two entries to be the first elements of the groups. Assign
            // each to a group.

            // debug code
            float initialArea = 0;
            
                Rectangle union = n.mbr.union(newRect);
                initialArea = union.area();
            

            System.Array.Copy(initialEntryStatus, 0, entryStatus, 0, maxNodeEntries);

            Node<T> newNode = null;
            newNode = new Node<T>(getNextNodeId(), n.level, maxNodeEntries);
            if (nodeMap.ContainsKey(newNode.nodeId) == false)
            {
                nodeMap.Add(newNode.nodeId, newNode);
            }
            else
            {
                nodeMap[newNode.nodeId] = newNode;
            }

            pickSeeds(n, newRect, newId, newNode); // this also sets the entryCount to 1

            // [Check if done] If all entries have been assigned, stop. If one
            // group has so few entries that all the rest must be assigned to it in 
            // order for it to have the minimum number m, assign them and stop. 
            while (n.entryCount + newNode.entryCount < maxNodeEntries + 1)
            {
                if (maxNodeEntries + 1 - newNode.entryCount == minNodeEntries)
                {
                    // assign all remaining entries to original openNode
                    for (int i = 0; i < maxNodeEntries; i++)
                    {
                        if (entryStatus[i] == ENTRY_STATUS_UNASSIGNED)
                        {
                            entryStatus[i] = ENTRY_STATUS_ASSIGNED;
                            n.mbr.add(n.entries[i]);
                            n.entryCount++;
                        }
                    }
                    break;
                }
                if (maxNodeEntries + 1 - n.entryCount == minNodeEntries)
                {
                    // assign all remaining entries to new openNode
                    for (int i = 0; i < maxNodeEntries; i++)
                    {
                        if (entryStatus[i] == ENTRY_STATUS_UNASSIGNED)
                        {
                            entryStatus[i] = ENTRY_STATUS_ASSIGNED;
                            newNode.addEntryNoCopy(n.entries[i], n.ids[i]);
                            n.entries[i] = null;
                        }
                    }
                    break;
                }

                // [Select entry to assign] Invoke algorithm pickNext to choose the
                // next entry to assign. Add it to the group whose covering rectangle 
                // will have to be enlarged least to accommodate it. Resolve ties
                // by adding the entry to the group with smaller area, then to the 
                // the one with fewer entries, then to either. Repeat from S2
                pickNext(n, newNode);
            }

            n.reorganize(this);

            // check that the MBR stored for each Node&lt;T&gt; is correct.
            if (INTERNAL_CONSISTENCY_CHECKING)
            {
                
            }

            // debug code
            

            return newNode;
        }

        /// <summary>
        /// Pick the seeds used to split a openNode.
        /// Select two entries to be the first elements of the groups
        /// </summary>
        /// <param name="n"></param>
        /// <param name="newRect"></param>
        /// <param name="newId"></param>
        /// <param name="newNode"></param>
        private void pickSeeds(Node<T> n, Rectangle newRect, int newId, Node<T> newNode)
        {
            // Find extreme rectangles along all dimension. Along each dimension,
            // find the entry whose rectangle has the highest low side, and the one 
            // with the lowest high side. Record the separation.
            float maxNormalizedSeparation = 0;
            int highestLowIndex = 0;
            int lowestHighIndex = 0;

            // for the purposes of picking seeds, take the MBR of the Node&lt;T&gt; to include
            // the new rectangle aswell.
            n.mbr.add(newRect);

        

            for (int d = 0; d < Rectangle.DIMENSIONS; d++)
            {
                float tempHighestLow = newRect.min[d];
                int tempHighestLowIndex = -1; // -1 indicates the new rectangle is the seed

                float tempLowestHigh = newRect.max[d];
                int tempLowestHighIndex = -1;

                for (int i = 0; i < n.entryCount; i++)
                {
                    float tempLow = n.entries[i].min[d];
                    if (tempLow >= tempHighestLow)
                    {
                        tempHighestLow = tempLow;
                        tempHighestLowIndex = i;
                    }
                    else
                    {  // ensure that the same index cannot be both lowestHigh and highestLow
                        float tempHigh = n.entries[i].max[d];
                        if (tempHigh <= tempLowestHigh)
                        {
                            tempLowestHigh = tempHigh;
                            tempLowestHighIndex = i;
                        }
                    }

                    // PS2 [Adjust for shape of the rectangle cluster] Normalize the separations
                    // by dividing by the widths of the entire set along the corresponding
                    // dimension
                    float normalizedSeparation = (tempHighestLow - tempLowestHigh) / (n.mbr.max[d] - n.mbr.min[d]);


                    // PS3 [Select the most extreme pair] Choose the pair with the greatest
                    // normalized separation along any dimension.
                    if (normalizedSeparation > maxNormalizedSeparation)
                    {
                        maxNormalizedSeparation = normalizedSeparation;
                        highestLowIndex = tempHighestLowIndex;
                        lowestHighIndex = tempLowestHighIndex;
                    }
                }
            }

            // highestLowIndex is the seed for the new openNode.
            if (highestLowIndex == -1)
            {
                newNode.addEntry(newRect, newId);
            }
            else
            {
                newNode.addEntryNoCopy(n.entries[highestLowIndex], n.ids[highestLowIndex]);
                n.entries[highestLowIndex] = null;

                // move the new rectangle into the space vacated by the seed for the new openNode
                n.entries[highestLowIndex] = newRect;
                n.ids[highestLowIndex] = newId;
            }

            // lowestHighIndex is the seed for the original openNode. 
            if (lowestHighIndex == -1)
            {
                lowestHighIndex = highestLowIndex;
            }

            entryStatus[lowestHighIndex] = ENTRY_STATUS_ASSIGNED;
            n.entryCount = 1;
            n.mbr.set(n.entries[lowestHighIndex].min, n.entries[lowestHighIndex].max);
        }



        
        /// <summary>
        /// Pick the next entry to be assigned to a group during a Node&lt;T&gt; split.
        /// [Determine cost of putting each entry in each group] For each 
        /// entry not yet in a group, calculate the area increase required
        /// in the covering rectangles of each group  
        /// </summary>
        /// <param name="n"></param>
        /// <param name="newNode"></param>
        /// <returns></returns>
        private int pickNext(Node<T> n, Node<T> newNode)
        {
            float maxDifference = float.NegativeInfinity;
            int next = 0;
            int nextGroup = 0;

            maxDifference = float.NegativeInfinity;


            for (int i = 0; i < maxNodeEntries; i++)
            {
                if (entryStatus[i] == ENTRY_STATUS_UNASSIGNED)
                {


                    float nIncrease = n.mbr.enlargement(n.entries[i]);
                    float newNodeIncrease = newNode.mbr.enlargement(n.entries[i]);
                    float difference = Math.Abs(nIncrease - newNodeIncrease);

                    if (difference > maxDifference)
                    {
                        next = i;

                        if (nIncrease < newNodeIncrease)
                        {
                            nextGroup = 0;
                        }
                        else if (newNodeIncrease < nIncrease)
                        {
                            nextGroup = 1;
                        }
                        else if (n.mbr.area() < newNode.mbr.area())
                        {
                            nextGroup = 0;
                        }
                        else if (newNode.mbr.area() < n.mbr.area())
                        {
                            nextGroup = 1;
                        }
                        else if (newNode.entryCount < maxNodeEntries / 2)
                        {
                            nextGroup = 0;
                        }
                        else
                        {
                            nextGroup = 1;
                        }
                        maxDifference = difference;
                    }
         
                }
            }

            entryStatus[next] = ENTRY_STATUS_ASSIGNED;

            if (nextGroup == 0)
            {
                n.mbr.add(n.entries[next]);
                n.entryCount++;
            }
            else
            {
                // move to new openNode.
                newNode.addEntryNoCopy(n.entries[next], n.ids[next]);
                n.entries[next] = null;
            }

            return next;
        }

        
        /// <summary>
        /// Recursively searches the tree for the nearest entry. Other queries
        /// call execute() on an IntProcedure when a matching entry is found; 
        /// however nearest() must store the entry Ids as it searches the tree,
        /// in case a nearer entry is found.
        /// Uses the member variable nearestIds to store the nearest
        /// entry IDs.
        /// </summary>
        /// <remarks>TODO rewrite this to be non-recursive?</remarks>
        /// <param name="p"></param>
        /// <param name="n"></param>
        /// <param name="nearestDistance"></param>
        /// <returns></returns>
        private float nearest(Point p, Node<T> n, float nearestDistance)
        {
            for (int i = 0; i < n.entryCount; i++)
            {
                float tempDistance = n.entries[i].distance(p);
                if (n.isLeaf())
                { // for leaves, the distance is an actual nearest distance 
                    if (tempDistance < nearestDistance)
                    {
                        nearestDistance = tempDistance;
                        nearestIds.Clear();
                    }
                    if (tempDistance <= nearestDistance)
                    {
                        nearestIds.Add(n.ids[i]);
                    }
                }
                else
                { // for index nodes, only go into them if they potentially could have
                    // a rectangle nearer than actualNearest
                    if (tempDistance <= nearestDistance)
                    {
                        // search the child openNode
                        nearestDistance = nearest(p, getNode(n.ids[i]), nearestDistance);
                    }
                }
            }
            return nearestDistance;
        }

    
        /// <summary>
        /// Recursively searches the tree for all intersecting entries.
        /// Immediately calls execute() on the passed IntProcedure when 
        /// a matching entry is found.
        /// [x] TODO rewrite this to be non-recursive? Make sure it
        /// doesn't slow it down.
        /// </summary>
        /// <param name="r"></param>
        /// <param name="v"></param>
        /// <param name="n"></param>
        private void intersects(Rectangle r, intproc v, Node<T> n)
        {
            for (int i = 0; i < n.entryCount; i++)
            {
                if (r.intersects(n.entries[i]))
                {
                    if (n.isLeaf())
                    {
                        v(n.ids[i]);
                    }
                    else
                    {
                        Node<T> childNode = getNode(n.ids[i]);
                        intersects(r, v, childNode);
                    }
                }
            }
        }

        /**
         * Used by delete(). Ensures that all nodes from the passed openNode
         * up to the root have the minimum number of entries.
         * 
         * Note that the parent and parentEntry stacks are expected to
         * contain the nodeIds of all parents up to the root.
         */

        private Rectangle oldRectangle = new Rectangle(0, 0, 0, 0);
        private IEnumerable<MGroup<T>> _getAllGroup;

        private void condenseTree(Node<T> l)
        {
            // CT1 [Initialize] Set n=l. Set the list of eliminated
            // nodes to be empty.
            Node<T> n = l;
            Node<T> parent = null;
            int parentEntry = 0;

            //TIntStack eliminatedNodeIds = new TIntStack();
            Stack<int> eliminatedNodeIds = new Stack<int>();

            // CT2 [Find parent entry] If N is the root, go to CT6. Otherwise 
            // let P be the parent of N, and let En be N's entry in P  
            while (n.level != treeHeight)
            {
                parent = getNode(parents.Pop());
                parentEntry = parentsEntry.Pop();

                // CT3 [Eliminiate under-full openNode] If N has too few entries,
                // delete En from P and add N to the list of eliminated nodes
                if (n.entryCount < minNodeEntries)
                {
                    parent.deleteEntry(parentEntry, minNodeEntries);
                    eliminatedNodeIds.Push(n.nodeId);
                }
                else
                {
                    // CT4 [Adjust covering rectangle] If N has not been eliminated,
                    // adjust EnI to tightly contain all entries in N
                    if (!n.mbr.Equals(parent.entries[parentEntry]))
                    {
                        oldRectangle.set(parent.entries[parentEntry].min, parent.entries[parentEntry].max);
                        parent.entries[parentEntry].set(n.mbr.min, n.mbr.max);
                        parent.recalculateMBR(oldRectangle);
                    }
                }
                // CT5 [Move up one level in tree] Set N=P and repeat from CT2
                n = parent;
            }

            // CT6 [Reinsert orphaned entries] Reinsert all entries of nodes in set Q.
            // Entries from eliminated leaf nodes are reinserted in tree leaves as in 
            // Insert(), but entries from higher level nodes must be placed higher in 
            // the tree, so that leaves of their dependent subtrees will be on the same
            // level as leaves of the main tree
            while (eliminatedNodeIds.Count > 0)
            {
                Node<T> e = getNode(eliminatedNodeIds.Pop());
                for (int j = 0; j < e.entryCount; j++)
                {
                    add(e.entries[j], e.ids[j], e.level);
                    e.entries[j] = null;
                }
                e.entryCount = 0;
                deletedNodeIds.Push(e.nodeId);
            }
        }

        /**
         *  Used by add(). Chooses a leaf to add the rectangle to.
         */
        private Node<T> chooseNode(Rectangle r, int level)
        {
            // CL1 [Initialize] Set N to be the root openNode
            Node<T> n = getNode(rootNodeId);
            parents.Clear();
            parentsEntry.Clear();

            // CL2 [Leaf check] If N is a leaf, return N
            while (true)
            {
        

                if (n.level == level)
                {
                    return n;
                }

                // CL3 [Choose subtree] If N is not at the desired level, let F be the entry in N 
                // whose rectangle FI needs least enlargement to include EI. Resolve
                // ties by choosing the entry with the rectangle of smaller area.
                float leastEnlargement = n.getEntry(0).enlargement(r);
                int index = 0; // index of rectangle in subtree
                for (int i = 1; i < n.entryCount; i++)
                {
                    Rectangle tempRectangle = n.getEntry(i);
                    float tempEnlargement = tempRectangle.enlargement(r);
                    if ((tempEnlargement < leastEnlargement) ||
                        ((tempEnlargement == leastEnlargement) &&
                         (tempRectangle.area() < n.getEntry(index).area())))
                    {
                        index = i;
                        leastEnlargement = tempEnlargement;
                    }
                }

                parents.Push(n.nodeId);
                parentsEntry.Push(index);

                // CL4 [Descend until a leaf is reached] Set N to be the child Node&lt;T&gt; 
                // pointed to by Fp and repeat from CL2
                n = getNode(n.ids[index]);
            }
        }

        /**
         * Ascend from a leaf Node&lt;T&gt; L to the root, adjusting covering rectangles and
         * propagating Node&lt;T&gt; splits as necessary.
         */
        private Node<T> adjustTree(Node<T> n, Node<T> nn)
        {
            // AT1 [Initialize] Set N=L. If L was split previously, set NN to be 
            // the resulting second openNode.

            // AT2 [Check if done] If N is the root, stop
            while (n.level != treeHeight)
            {

                // AT3 [Adjust covering rectangle in parent entry] Let P be the parent 
                // Node<T> of N, and let En be N's entry in P. Adjust EnI so that it tightly
                // encloses all entry rectangles in N.
                Node<T> parent = getNode(parents.Pop());
                int entry = parentsEntry.Pop();

     
                if (!parent.entries[entry].Equals(n.mbr))
                {
                    parent.entries[entry].set(n.mbr.min, n.mbr.max);
                    parent.mbr.set(parent.entries[0].min, parent.entries[0].max);
                    for (int i = 1; i < parent.entryCount; i++)
                    {
                        parent.mbr.add(parent.entries[i]);
                    }
                }

                // AT4 [Propagate Node<T> split upward] If N has a partner NN resulting from 
                // an earlier split, create a new entry Enn with Ennp pointing to NN and 
                // Enni enclosing all rectangles in NN. Add Enn to P if there is room. 
                // Otherwise, invoke splitNode to produce P and PP containing Enn and
                // all P's old entries.
                Node<T> newNode = null;
                if (nn != null)
                {
                    if (parent.entryCount < maxNodeEntries)
                    {
                        parent.addEntry(nn.mbr, nn.nodeId);
                    }
                    else
                    {
                        newNode = splitNode(parent, nn.mbr.copy(), nn.nodeId);
                    }
                }

                // AT5 [Move up to next level] Set N = P and set NN = PP if a split 
                // occurred. Repeat from AT2
                n = parent;
                nn = newNode;

                parent = null;
                newNode = null;
            }

            return nn;
        }

        /**
         * Check the consistency of the tree.
         */
        private void checkConsistency(int nodeId, int expectedLevel, Rectangle expectedMBR)
        {
            // go through the tree, and check that the internal data structures of 
            // the tree are not corrupted.    
            Node<T> n = getNode(nodeId);

       

            Rectangle calculatedMBR = calculateMBR(n);

    
            // Check for corruption where a parent entry is the same object as the child MBR
      

            for (int i = 0; i < n.entryCount; i++)
            {

                if (n.level > 1)
                { // if not a leaf
                    checkConsistency(n.ids[i], n.level - 1, n.entries[i]);
                }
            }
        }

        /**
         * Given a Node<T> object, calculate the Node<T> MBR from it's entries.
         * Used in consistency checking
         */
        private Rectangle calculateMBR(Node<T> n)
        {
            Rectangle mbr = new Rectangle(n.entries[0].min, n.entries[0].max);

            for (int i = 1; i < n.entryCount; i++)
            {
                mbr.add(n.entries[i]);
            }
            return mbr;
        }

        public int Count
        {
            get
            {
                return this.msize;
            }
        }
    }
}