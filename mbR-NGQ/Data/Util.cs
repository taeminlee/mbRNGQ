using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mbR_NGQ.Data
{
    public static class Util
    {
        public static IEnumerable<IEnumerable<T>> Combinations<T>(this IEnumerable<T> elements, int k)
        {
            // 0,1,2 => 0,1 => 0,2 => 1,2 => 0 => 1 => 2
            return k == 0 ? new[] { new T[0] } :
              elements.SelectMany((e, i) =>
                Combinations(elements.Skip(i + 1), k - 1).Select(c => (new[] { e }).Concat(c)));
        }
    }
}
