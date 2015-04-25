using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace mbR_NGQ.UI
{
    public partial class DebugForm : Form
    {
        public DebugForm()
        {
            InitializeComponent();
        }

        public void DebugWrite(string str)
        {
            richTextBox1.AppendText(str + Environment.NewLine);
        }
    }
}
