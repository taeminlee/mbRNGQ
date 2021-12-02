using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Runtime.InteropServices.ComTypes;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace mbR_NGQ
{
    public partial class OptionForm : Form
    {
        public OptionForm()
        {
            InitializeComponent();
        }

        private void textBox1_TextChanged(object sender, EventArgs e)
        {
            int m;
            if (int.TryParse(this.textBox1.Text, out m) == false)
            {
                MessageBox.Show("Incorrect Value", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                OptionForm_Load(this, null);
            }
            else
            {
                Config.m = m;
            }
        }

        private void textBox2_TextChanged(object sender, EventArgs e)
        {
            int k;
            if (int.TryParse(this.textBox2.Text, out k) == false)
            {
                MessageBox.Show("Incorrect Value", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                OptionForm_Load(this, null);
            }
            else
            {
                Config.k = k;
            }
        }

        private void textBox3_TextChanged(object sender, EventArgs e)
        {
            int n;
            if (int.TryParse(this.textBox3.Text, out n) == false)
            {
                MessageBox.Show("Incorrect Value", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                OptionForm_Load(this, null);
            }
            else
            {
                Config.n = n;
            }
        }

        private void textBox4_TextChanged(object sender, EventArgs e)
        {
            float minX;
            if (float.TryParse(this.textBox4.Text, out minX) == false)
            {
                MessageBox.Show("Incorrect Value", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                OptionForm_Load(this, null);
            }
            else
            {
                Config.minX = minX;
            }
        }

        private void textBox5_TextChanged(object sender, EventArgs e)
        {
            float maxX;
            if (float.TryParse(this.textBox5.Text, out maxX) == false)
            {
                MessageBox.Show("Incorrect Value", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                OptionForm_Load(this, null);
            }
            else
            {
                Config.maxX = maxX;
            }
        }

        private void textBox6_TextChanged(object sender, EventArgs e)
        {
            float minY;
            if (float.TryParse(this.textBox6.Text, out minY) == false)
            {
                MessageBox.Show("Incorrect Value", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                OptionForm_Load(this, null);
            }
            else
            {
                Config.minY = minY;
            }
        }

        private void textBox7_TextChanged(object sender, EventArgs e)
        {
            float maxY;
            if (float.TryParse(this.textBox7.Text, out maxY) == false)
            {
                MessageBox.Show("Incorrect Value", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                OptionForm_Load(this, null);
            }
            else
            {
                Config.maxY = maxY;
            }
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            //Config.DataGenerateMethods = (DataSetGenerateMethods)Enum.Parse(typeof (DataSetGenerateMethods), comboBox1.Text);
        }

        private void OptionForm_Load(object sender, EventArgs e)
        {
            this.textBox1.Text = Config.m.ToString();
            this.textBox2.Text = Config.k.ToString();
            this.textBox3.Text = Config.n.ToString();
            this.textBox4.Text = Config.minX.ToString();
            this.textBox5.Text = Config.maxX.ToString();
            this.textBox6.Text = Config.minY.ToString();
            this.textBox7.Text = Config.maxY.ToString();
            /*foreach (DataSetGenerateMethods method in Enum.GetValues(typeof(DataSetGenerateMethods)))
            {
                int id = comboBox1.Items.Add(method.ToString());
            }
            this.comboBox1.SelectedItem = this.comboBox1.Items[(int) Config.DataGenerateMethods];*/
        }

        private void button1_Click(object sender, EventArgs e)
        {
            RTreeViewer rTreeViewer = new RTreeViewer();
            rTreeViewer.Show();
        }
    }
}
