using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace GeneticsProject
{
    public struct GeneticData
    {
        public string name; 
        public string organism;
        public string formula; 
    }

    class Program
    {
        static List<GeneticData> data = new List<GeneticData>();
        static int count = 1;
        const bool behaivour_okey = true;
        const bool behaivour_no_okey = false;
        static string GetFormula(string proteinName)
        {
            foreach (GeneticData item in data)
            {
                if (item.name.Equals(proteinName)) return item.formula;
            }
            return null;
        }
        static void ReadGeneticData(string filename)
        {
            StreamReader reader = new StreamReader(filename);
            while (!reader.EndOfStream)
            {
                string line = reader.ReadLine();
                string[] fragments = line.Split('\t');
                GeneticData protein;
                protein.name = fragments[0];
                protein.organism = fragments[1];
                protein.formula = fragments[2];
                if (IsValid(Decoding(protein.formula)) == behaivour_no_okey)
                {
                    throw new Exception($"The formula <{protein.name}> was not volidated");
                }
                data.Add(protein);
                count++;
            }
            reader.Close();
        }
        static void ReadHandleCommands(string filename, StreamWriter writer)
        {
            StreamReader reader = new StreamReader(filename);
            int counter = 0;
            while (!reader.EndOfStream)
            {
                string line = reader.ReadLine(); counter++;
                string[] command = line.Split('\t');
                string doing = command[0];
                switch (doing)
                {
                    case "search":
                        writer.WriteLine($"{counter.ToString("D3")}\t{"search"}\t{Decoding(command[1])}");
                        int index = Search(command[1]);
                        if (index != -1)
                        {
                            writer.WriteLine("organism\t\t\tprotein");
                            writer.WriteLine($"{data[index].organism}\t\t{data[index].name}");
                            writer.WriteLine("------------------------------------------------------------------------------");
                        }
                        else
                        {
                            writer.WriteLine("organism\t\t\t\tprotein");
                            writer.WriteLine("NOT FOUND");
                            writer.WriteLine("------------------------------------------------------------------------------");
                        }
                        break;
                    case "diff":
                        writer.WriteLine($"{counter.ToString("D3")}\tdiff\t{command[1]}\t{command[2]}");
                        writer.WriteLine("amino-acids difference:");
                        if (Diff(command[1], command[2]) >= 0)
                        {
                            writer.WriteLine(Diff(command[1], command[2]));
                            writer.WriteLine("------------------------------------------------------------------------------");
                        }
                        else if (Diff(command[1], command[2]) == -1)
                        {
                            writer.WriteLine($"MISSING : {command[1]}");
                            writer.WriteLine("------------------------------------------------------------------------------");
                        }
                        else if (Diff(command[1], command[2]) == -2)
                        {
                            writer.WriteLine($"MISSING : {command[2]}");
                            writer.WriteLine("------------------------------------------------------------------------------");
                        }
                        break;
                    case "mode":
                        writer.WriteLine($"{counter.ToString("D3")}\t{"mode"}\t{command[1]}");
                        writer.WriteLine("amino-acid occurs: ");
                        Mode(command[1], out string letter, out int a);
                        if (a == 0)
                        {
                            writer.WriteLine($"MISSING : {command[1]}");
                            writer.WriteLine("------------------------------------------------------------------------------");
                        }
                        else
                        {
                            writer.WriteLine($"{letter}\t\t{a}");
                            writer.WriteLine("------------------------------------------------------------------------------");
                        }
                        break;
                    default:
                        throw new Exception("Wrong command");
                }
            }
            reader.Close();
            writer.Close();
        }
        static bool IsValid(string formula)
        {
            List<char> letters = new List<char>() { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'B' };
            foreach (char ch in formula)
            {
                if (!letters.Contains(ch)) return false;
            }
            return true;
        }
        static string ENcoding(string formula)
        {
            string encoded = String.Empty;
            for (int i = 0; i < formula.Length; i++)
            {
                char ch = formula[i];
                int count = 1;
                while (i < formula.Length - 1 && formula[i + 1] == ch)
                {
                    count++;
                    i++;
                }
                if (count > 2 && count <= 9) encoded = encoded + count + ch;
                if (count == 1) encoded = encoded + ch;
                if (count == 2) encoded = encoded + ch + ch;
                if (count > 9) throw new Exception("More than 9 consecutive letter");
            }
            return encoded;
        }
        static string Decoding(string formula)
        {
            string decoded = String.Empty;
            for (int i = 0; i < formula.Length; i++)
            {
                if (char.IsDigit(formula[i]))
                {
                    char letter = formula[i + 1];
                    int conversion = formula[i] - '0';
                    for (int j = 0; j < conversion - 1; j++) decoded = decoded + letter;
                }
                else decoded = decoded + formula[i];
            }
            return decoded;
        }
        static int Search(string amino_acid)
        {
            string decoded = Decoding(amino_acid);
            for (int i = 0; i < data.Count; i++)
            {
                if (Decoding(data[i].formula).Contains(decoded)) return i;
            }
            return -1;
        }
        static int Diff(string protein1, string protein2)
        {

            string acid1 = string.Empty;
            bool behaviour_first = false; bool behaivour_second = false;
            string acid2 = string.Empty;
            int max_of_length = 0;
            int diff_of_length = 0;
            for (int i = 0; i< data.Count; i++)
            {
                if (data[i].name == protein1) 
                {
                    acid1 = Decoding(GetFormula(protein1));
                    behaviour_first = behaivour_okey;
                }
                if (data[i].name == protein2) 
                {
                    acid2 = Decoding(GetFormula(protein2));
                    behaivour_second = behaivour_okey; 
                }
            }
            if (acid1.Length >= acid2.Length) 
            {
                max_of_length =acid1.Length; 
                diff_of_length = max_of_length - acid2.Length;
            } 
            else 
            { 
                max_of_length = acid2.Length; diff_of_length = max_of_length - acid1.Length; 
            }
            int number = 0 + diff_of_length; 
            for (int i = 0; i < max_of_length - diff_of_length; i++) 
            {
                if (acid1[i] != acid2[i]) number++;

            }
            if (behaviour_first == behaivour_okey && behaivour_second == behaivour_okey) 
            { 
                return number;
            }
            else if (behaviour_first == behaivour_no_okey)
            { 
                return -1;
            }
            else
            {
                return -2;
            }
        }
        static void Mode(string protein, out string letter, out int a) 
        {
            string formula = GetFormula(protein);
            if (formula == null)
            {
                 a = 0;
                letter = null;
            }
            else
            {
                Dictionary<string, int> dict = new Dictionary<string, int>();
                char[] check = (Decoding(formula)).ToCharArray();
                for (int i = 0; i < check.Length; i++)
                {
                    if (dict.ContainsKey(Convert.ToString(check[i])))
                    {
                        dict[Convert.ToString(check[i])]++;
                    }
                    else
                    {
                        dict.Add(Convert.ToString(check[i]), 1);
                    }
                }
                int[] keys = dict.Values.ToArray();
                Array.Sort(keys);
                int m = keys[keys.Length - 1];
                var max_num = new KeyValuePair<string, int>("Z", m);
                foreach (var item in dict)
                {
                    if (item.Value == max_num.Value)
                    {
                        if (Convert.ToChar(item.Key) < Convert.ToChar(max_num.Key))
                        {
                            max_num = item;
                        }
                    }

                }
                letter = max_num.Key;
                a = Convert.ToInt32(max_num.Value);
            }
        }
        
        static void Main(string[] args)
        {
            try
            {
                File.Delete("output.txt");
                ReadGeneticData("sequences.0.txt");
                Console.WriteLine(ENcoding("AAAAAAAAATATTTCGCTTTTCAAAAATTGTCAGATGAGAGAAAAAATAAAA"));
                FileStream f = new FileStream("OutPut.txt", FileMode.OpenOrCreate);
                StreamWriter writer = new StreamWriter(f, Encoding.Default);
                writer.WriteLine("Lisovskyi Kirill\nGenetic Searching");
                writer.WriteLine("------------------------------------------------------------------------------");
                ReadHandleCommands("commands.0.txt",writer);
        }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
}
    }
}
