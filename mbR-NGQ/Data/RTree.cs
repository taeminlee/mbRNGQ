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
using System.Drawing.Drawing2D;
using System.Linq;
using System.Net;
using System.Security.Permissions;
using System.Text.RegularExpressions;
using mbR_NGQ;
using mbR_NGQ.Data;
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

        // to visualize mbr easily... , 20150105 ���¹�
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
        private volatile int idcounter = 0; // to know what is assigned easily... 20150106 ���¹�

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
        ///in a node. The default value is 10, which is used if the property is
        ///not specified, or is less than 2.</param>
        /// <param name="MinNodeEntries">This specifies the minimum number of entries
        ///in a node. The default value is half of the MaxNodeEntries value (rounded
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
            // per node, but will be inefficient.
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

            // I2 [Add record to leaf node] If L has room for another entry, 
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
                nodeMap.Add(rootNodeId, root);
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
            // entry is not a leaf node, delete the root (it's entry becomes the new root)
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

            public Rectangle GetCombinationMBR(double delta, Point queryPoint)
            {
                
                Rectangle mbr = new Rectangle(Config.minX, Config.minY, Config.maxX, Config.maxY);
                mbr = UnionLocalDeltaMBR(mbr, queryPoint, delta);

                foreach (Object obj in this._elems)
                {
                    if (obj == null)
                        continue;
                    else
                    {
                        if (obj.GetType() == typeof (Point))
                        {
                            mbr = UnionLocalDeltaMBR(mbr, (Point)obj, delta);
                        }
                        else if (obj.GetType() == typeof (Node<T>))
                        {
                            mbr = UnionLocalDeltaMBR(mbr, (Node<T>)obj, delta);
                        }
                    }
                }
                
                return mbr;
            }

            private Rectangle UnionLocalDeltaMBR(Rectangle mbr, Point point, double delta)
            {
                mbr.min[0] = Math.Max(mbr.min[0], point.coordinates[0] - (float)delta);
                mbr.min[1] = Math.Max(mbr.min[1], point.coordinates[1] - (float)delta);
                mbr.max[0] = Math.Min(mbr.max[0], point.coordinates[0] + (float)delta);
                mbr.max[1] = Math.Min(mbr.max[1], point.coordinates[1] + (float)delta);

                return mbr;
            }

            private Rectangle UnionLocalDeltaMBR(Rectangle mbr, Node<T> node, double delta)
            {
                mbr.min[0] = Math.Max(mbr.min[0], node.mbr.min[0] - (float)delta);
                mbr.min[1] = Math.Max(mbr.min[1], node.mbr.min[1] - (float)delta);
                mbr.max[0] = Math.Min(mbr.max[0], node.mbr.max[0] + (float)delta);
                mbr.max[1] = Math.Min(mbr.max[1], node.mbr.max[1] + (float)delta);

                return mbr;
            }

            public MGroup(object[] _elems)
            {
                System.Diagnostics.Debug.Assert(_elems.Length == Config.m);
                this._elems = _elems;
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
                    double interDist = GetMaxMinInterDist(r, q);
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
                    for (int j = i + 1; j < _elems.Length; j++)
                    {
                        if (_elems[i] == _elems[j]) continue;

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
                return newMGroup;
            }
        }

        private double globalMaxDist;
        private double delta;
        private int KCounter = 0;
        private List<MGroup<T>> candidateGroups;
        private Dictionary<object, List<MGroup<T>>> invertedCandidateGroups;
        private MGroup<T>[] resultGroups;
        private RTree<object>[] invertedElemRTrees;
        private Dictionary<object, Rectangle> elemRecDic;
        private double[] thetaPool;
        private MGroup<T>[] thetaGroups;
        private List<object> objPool;
        private List<object>[] invertedObj;
        private int leafCount;
        private BitArray groupBitArray;
        private int accessCount;

        public List<MGroup<T>> NearestGroup(Point q)
        {
            List<MGroup<T>> retval = new List<MGroup<T>>();

            Node<T> rootNode = getNode(rootNodeId);

            globalMaxDist = new Point(Config.minX, Config.minY, 0).GetDist(new Point(Config.maxX, Config.maxY, 0));
            delta = globalMaxDist;
            thetaPool = new double[Config.k];
            thetaGroups = new MGroup<T>[Config.k];
            candidateGroups = new List<MGroup<T>>();
            invertedCandidateGroups = new Dictionary<object, List<MGroup<T>>>();
            elemRecDic = new Dictionary<object, Rectangle>();
            resultGroups = new MGroup<T>[Config.k];
            objPool = new List<object>();
            invertedElemRTrees = new RTree<object>[Config.m];
            invertedObj = new List<object>[Config.m];
            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                invertedObj[mIdx] = new List<object>();
                invertedElemRTrees[mIdx] = new RTree<object>();
            }
            KCounter = 0;
            for (int i = 0; i < thetaPool.Length; i++)
            {
                thetaPool[i] = globalMaxDist;
            }
            invertedObj = new List<object>[Config.m];
            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                invertedObj[mIdx] = new List<object>();
            }
            leafCount = 0;
            accessCount = 0;
            groupBitArray = new BitArray(Config.m);
            groupBitArray.SetAll(true);

            bool allSet = true;
            MGroup<T> rootGroup = new MGroup<T>(new object[Config.m]);
            Rectangle rootRec = new Rectangle(rootNode as Node<Point>);
            for (int idx = 0; idx < Config.m; idx++)
            {
                rootGroup.Elems[idx] = rootNode;
                invertedElemRTrees[idx].Add(rootRec, rootNode);
            }
            System.Diagnostics.Debug.Assert(allSet == true);
            candidateGroups.Add(rootGroup);
            invertedCandidateGroups[rootNode] = new List<MGroup<T>>();
            invertedCandidateGroups[rootNode].Add(rootGroup);
            NearestGroup(q, rootNode);

            var v = resultGroups.ToList();
            v.Sort();
            resultGroups = v.ToArray();

            foreach (MGroup<T> ng in resultGroups)
            {
                retval.Add(ng);
            }

            return retval;
        }

        private void NearestGroup(Point q, Node<T> node)
        {
            accessCount++;
            node.visited = true;

            Rectangle rec = new Rectangle(node as Node<Point>);
            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                if (node.bitArray[mIdx] == true)
                {
                    invertedElemRTrees[mIdx].Delete(rec, node);
                }
            }

            //1. Initialize Inverted List
            //       [NOT NOW]      which Sorted order by coordinate

            List<object> newObjPool = new List<object>();
            List<object>[] newInvertedObj = new List<object>[Config.m];
            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                newInvertedObj[mIdx] = new List<object>();
            }

            InsertOpenElements(q, node, newInvertedObj, newObjPool);

            // Make Combination Group (Computation Cost)

            MGroup<T>[] openGroups = invertedCandidateGroups[node].ToArray();

            foreach (MGroup<T> curGroup in openGroups)
            {
                candidateGroups.Remove(curGroup);
                MakeCombinationGroup(q, curGroup, node, invertedElemRTrees, newObjPool, newInvertedObj);
            }
            invertedCandidateGroups.Remove(node);

            candidateGroups.Sort();

            RemoveGroupByTheta(q);

            int currentKCounter = KCounter;

            Node<T> bestNode = null;
            double bestDist = globalMaxDist;

            List<MGroup<T>> resultAddedGroups = new List<MGroup<T>>();

            for (int idx = 0; idx < Math.Min(candidateGroups.Count, Config.k - currentKCounter); idx++)
            {
                MGroup<T> group = candidateGroups[idx];
                if (bestNode == null && group.IsLeaf())
                {
                    resultGroups[KCounter++] = group;

                    if (KCounter == Config.k)
                        return;

                    DeQueueDeltaPool();
                    resultAddedGroups.Add(group);
                }
                else
                {
                    if (bestNode != null) continue;
                    foreach (object obj in group.Elems)
                    {
                        if (obj.GetType() == typeof (Node<T>))
                        {
                            Node<T> innerNode = (Node<T>)obj;
                            if (innerNode.mbr.distance(q) < bestDist)
                            {
                                bestDist = innerNode.mbr.distance(q);
                                bestNode = @innerNode;
                            }
                        }
                    }
                    break;
                }
            }

            foreach (MGroup<T> group in resultAddedGroups)
            {
                candidateGroups.Remove(group);
            }

            System.Diagnostics.Debug.Assert(bestNode != null);
            System.Diagnostics.Debug.Assert(invertedCandidateGroups.ContainsKey(bestNode) != false);
            
            NearestGroup(q, bestNode);
        }

        private void DeQueueDeltaPool()
        {
            List<double> thetaPoolList = thetaPool.ToList();
            List<MGroup<T>> thetaGroupList = thetaGroups.ToList();
            int thetaIdx = thetaPool.ToList().IndexOf(thetaPool.Min());
            thetaPoolList.RemoveAt(thetaIdx);
            thetaGroupList.RemoveAt(thetaIdx);
            thetaPool = thetaPoolList.ToArray();
            thetaGroups = thetaGroupList.ToArray();
        }

        private void InsertOpenElements(Point q, Node<T> node, List<object>[] newInvertedObj, List<object> newObjPool)
        {
            // Insert child nodes into inverted List
            if (node.isLeaf() == false)
            {
                Node<T>[] nodes = new Node<T>[node.entryCount];
                for (int i = 0; i < node.entryCount; i++)
                {
                    nodes[i] = NodeMap[node.ids[i]];
                }
                foreach (Node<T> childNode in nodes)
                {
                    Rectangle elemRec = new Rectangle(childNode as Node<Point>);
                    elemRecDic.Add(childNode, elemRec);

                    if (childNode.mbr.distance(q) > delta) continue;
                    bool allSet = true;
                    for (int mIdx = 0; mIdx < Config.m; mIdx++)
                    {
                        if (childNode.bitArray[mIdx] == true)
                        {
                            invertedElemRTrees[mIdx].Add(elemRec, childNode);
                            newInvertedObj[mIdx].Add(childNode);
                        }
                        else
                            allSet = false;
                    }
                    newObjPool.Add(childNode);
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
                    System.Diagnostics.Debug.Assert(obj.GetType() == typeof (Point));
                    if (obj.GetType() != typeof (Point)) continue;

                    Point point = obj as Point;

                    if (point.GetDist(q) > delta) continue;

                    Rectangle elemRec = new Rectangle(point);
                    elemRecDic.Add(point, elemRec);
                    invertedElemRTrees[point.category].Add(elemRec, point);

                    newInvertedObj[point.category].Add(point);
                    newObjPool.Add(point);
                }
            }
        }

        private void MakeCombinationGroup(Point queryPoint, MGroup<T> openGroup, Node<T> openNode, RTree<object>[] invertedElemRTree, List<object> newElemPool, List<object>[] newInvertedElem)
        {
            if(openGroup.IsDictator() && !openNode.isLeaf())
            { 
                MakeCombinationGroup(queryPoint, newElemPool, newInvertedElem);
                return;
            }

            MGroup<T> templateGroup = new MGroup<T>(new object[Config.m]);
            if(openGroup.IsDictator() == false)
            { 
                for (int mIdx = 0; mIdx < Config.m; mIdx++)
                {
                    if (openGroup.Elems[mIdx] != openNode)
                        templateGroup.Elems[mIdx] = openGroup.Elems[mIdx];
                }
            }

            MakeCombinationGroup(queryPoint, templateGroup, invertedElemRTree);
        }

        // When the OpenGroup is not a Dictator Group, Make Combination Groups using delta constraint rule
        private void MakeCombinationGroup(Point queryPoint, MGroup<T> templateGroup, RTree<object>[] invertedElemRTree)
        {
            List<MGroup<T>> partialGroups = new List<MGroup<T>>();

            partialGroups.Add(templateGroup);

            while (partialGroups.Count != 0)
            {
                List<MGroup<T>> removePartialGroups = new List<MGroup<T>>();
                List<MGroup<T>> addPartialGroups = new List<MGroup<T>>();

                foreach (MGroup<T> partialGroup in partialGroups)
                {
                    if (partialGroup.IsAllSet())
                    {
                        AddGroupToCandidate(queryPoint, partialGroup);
                        removePartialGroups.Add(partialGroup);
                        continue;
                    }

                    Rectangle searchArea = partialGroup.GetCombinationMBR(delta, queryPoint);

                    int idx = partialGroup.GetFirstNullIndex();

                    System.Diagnostics.Debug.Assert(idx != -1);

                    foreach (object element in invertedElemRTree[idx].Contains(searchArea))
                    {
                        System.Diagnostics.Debug.Assert(element != null);

                        MGroup<T> newGroup = partialGroup.Copy();
                        AddElementToGroup(element, newGroup);
                        addPartialGroups.Add(newGroup);
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
        }

        // When the OpenGroup is a Dictator Group and openNode is internal Node, Make Combination Groups using newly added objets only.
        private void MakeCombinationGroup(Point queryPoint, List<object> newElemPool, List<object>[] newInvertedElem )
        {
            List<MGroup<T>> partialGroups = new List<MGroup<T>>();

            foreach (object elem in newElemPool)
            {
                partialGroups.Add(GetNewGroupByElement(elem));
            }

            List < MGroup < T >> removeGroups = new List<MGroup<T>>();

            foreach (MGroup<T> group in partialGroups)
            {
                if (group.IsAllSet())
                {
                    AddGroupToCandidate(queryPoint, group);
                    removeGroups.Add(group);
                }
            }

            foreach (MGroup<T> group in removeGroups)
            {
                partialGroups.Remove(group);
            }

            while(partialGroups.Count != 0)
            {
                List<MGroup<T>> removePartialGroups = new List<MGroup<T>>();
                List<MGroup<T>> addPartialGroups = new List<MGroup<T>>();

                foreach (MGroup<T> partialGroup in partialGroups)
                {
                    if (partialGroup.IsAllSet())
                    {
                        AddGroupToCandidate(queryPoint, partialGroup);
                        removePartialGroups.Add(partialGroup);
                        continue;
                    }

                    int idx = partialGroup.GetFirstNullIndex();

                    System.Diagnostics.Debug.Assert(idx != -1);

                    foreach (object element in newInvertedElem[idx])
                    {
                        MGroup<T> newGroup = partialGroup.Copy();
                        AddElementToGroup(element, newGroup);
                        addPartialGroups.Add(newGroup);
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
        }

        private bool AddGroupToCandidate(Point queryPoint, MGroup<T> partialGroup)
        {
            if (DuplicateCheck(partialGroup) == false) return false;

            partialGroup.GetGroupMinDist(queryPoint);
            if (partialGroup.MinDist > delta)
                return false;
            partialGroup.GetGroupMaxDist(queryPoint);
            UpdateDelta(partialGroup);
            candidateGroups.Add(partialGroup);
            return true;
        }

        private bool DuplicateCheck(MGroup<T> partialGroup)
        {
            foreach (MGroup<T> group in candidateGroups)
            {
                bool checkBit = false;
                for (int mIdx = 0; mIdx < Config.m; mIdx++)
                {
                    if (partialGroup.Elems[mIdx] != group.Elems[mIdx])
                    {
                        checkBit = true;
                    }
                }
                if (checkBit == false)
                    return checkBit;
            }
            return true;
        }

        private MGroup<T> GetNewGroupByElement(object element)
        {
            MGroup<T> group = new MGroup<T>(new object[Config.m]);
            AddElementToGroup(element, @group);
            return group;
        }

        private void AddElementToGroup(object element, MGroup<T> @group)
        {
            if (element.GetType() == typeof (Node<T>))
            {
                Node<T> node = element as Node<T>;

                if(invertedCandidateGroups.ContainsKey(node) == false)
                    invertedCandidateGroups[node] = new List<MGroup<T>>();
                invertedCandidateGroups[node].Add(@group);

                for (int mIdx = 0; mIdx < Config.m; mIdx++)
                {
                    if (node.bitArray[mIdx] == true)
                        @group.Elems[mIdx] = node;
                }
            }
            else if (element.GetType() == typeof (Point))
            {
                Point obj = element as Point;
                @group.Elems[obj.category] = obj;
            }
        }

        private int getIntFromBitArray(BitArray bitArray)
        {

            if (bitArray.Length > 32)
                throw new ArgumentException("Argument length shall be at most 32 bits.");

            int[] array = new int[1];
            bitArray.CopyTo(array, 0);
            return array[0];

        }

        private void RemoveGroupByTheta(Point q)
        {
            List<MGroup<T>> removeGroups = new List<MGroup<T>>();
            List<object> removeKeys = new List<object>();
            foreach (var group in candidateGroups)
            {
                if (group.MinDist > delta)
                {
                    removeGroups.Add(group);
                }
            }
            foreach (var group in removeGroups)
            {
                candidateGroups.Remove(group);

                foreach (KeyValuePair<object, List<MGroup<T>>> kv in invertedCandidateGroups)
                {
                    if (kv.Key.GetType() == typeof(Node<T>))
                    {
                        Node<T> innerNode = kv.Key as Node<T>;
                        if (innerNode.mbr.distance(q) > delta)
                        {
                            kv.Value.Clear();
                            removeKeys.Add(kv.Key);
                            continue;
                        }
                    }
                    else if (kv.Key.GetType() == typeof(Point))
                    {
                        Point innerPoint = kv.Key as Point;
                        if (innerPoint.GetDist(q) > delta)
                        {
                            kv.Value.Clear();
                            removeKeys.Add(kv.Key);
                            continue;
                        }
                    }
                    if (kv.Value.Contains(group))
                    {
                        kv.Value.Remove(group);
                    }
                    if(kv.Value.Count == 0)
                        removeKeys.Add(kv.Key);
                }
            }
            foreach (object obj in removeKeys)
            {
                invertedCandidateGroups.Remove(obj);
            }
        }

        private void UpdateInvertedCandidateGroups(MGroup<T> newGroup)
        {
            List<object> skipList = new List<object>();
            foreach (object obj in newGroup.Elems)
            {
                if (skipList.Contains(obj)) continue;
                skipList.Add(obj);
                if (obj.GetType() == typeof (Node<T>) == false) continue;
                if (invertedCandidateGroups.ContainsKey(obj) == false)
                    invertedCandidateGroups[obj] = new List<MGroup<T>>();
                invertedCandidateGroups[obj].Add(newGroup);
            }
        }

        private void UpdateDelta(MGroup<T> newGroup)
        {
            if (thetaPool.Max() > newGroup.MaxDist)
            {
                int idx = thetaPool.ToList().IndexOf(thetaPool.Max());
                thetaPool[idx] = newGroup.MaxDist;
                delta = thetaPool.Max();
                thetaGroups[idx] = newGroup;
            }
        }

        private IEnumerable<MGroup<T>> MakeCombinationGroup(MGroup<T> openGroup, Node<T> openNode)
        {
            List<MGroup<T>> groups = new List<MGroup<T>>();
            List<MGroup<T>> partialGroups = new List<MGroup<T>>();

            MGroup<T> templateGroup = new MGroup<T>(new object[Config.m]);
            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                if(openGroup.Elems[mIdx] != openNode)
                    templateGroup.Elems[mIdx] = openGroup.Elems[mIdx];
            }

            partialGroups.Add(templateGroup);

            for (int mIdx = 0; mIdx < Config.m; mIdx++)
            {
                if (templateGroup.Elems[mIdx] != null) continue;

                List<MGroup<T>> removePartialGroups = new List<MGroup<T>>();
                List<MGroup<T>> addPartialGroups = new List<MGroup<T>>();
                foreach (object obj in invertedObj[mIdx])
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
                        groups.Add(group);
                    else
                        partialGroups.Add(group);
                }
            }

            return groups;
        }

        private List<MGroup<T>> goldGroups;


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
                    this.NearestGroup(queryPoint);
                    return false;
                }
            }

            return true;
        }


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
            // MBR of the root node. If no intersection, return immediately.

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
        /// Get a Node&lt;T&gt; object, given the ID of the node.
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
        /// Split a node. Algorithm is taken pretty much verbatim from
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
            nodeMap.Add(newNode.nodeId, newNode);

            pickSeeds(n, newRect, newId, newNode); // this also sets the entryCount to 1

            // [Check if done] If all entries have been assigned, stop. If one
            // group has so few entries that all the rest must be assigned to it in 
            // order for it to have the minimum number m, assign them and stop. 
            while (n.entryCount + newNode.entryCount < maxNodeEntries + 1)
            {
                if (maxNodeEntries + 1 - newNode.entryCount == minNodeEntries)
                {
                    // assign all remaining entries to original node
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
                    // assign all remaining entries to new node
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
        /// Pick the seeds used to split a node.
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

            // highestLowIndex is the seed for the new node.
            if (highestLowIndex == -1)
            {
                newNode.addEntry(newRect, newId);
            }
            else
            {
                newNode.addEntryNoCopy(n.entries[highestLowIndex], n.ids[highestLowIndex]);
                n.entries[highestLowIndex] = null;

                // move the new rectangle into the space vacated by the seed for the new node
                n.entries[highestLowIndex] = newRect;
                n.ids[highestLowIndex] = newId;
            }

            // lowestHighIndex is the seed for the original node. 
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
                // move to new node.
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
                        // search the child node
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
         * Used by delete(). Ensures that all nodes from the passed node
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

                // CT3 [Eliminiate under-full node] If N has too few entries,
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
            // CL1 [Initialize] Set N to be the root node
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
            // the resulting second node.

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