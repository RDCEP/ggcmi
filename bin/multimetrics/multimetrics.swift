type file;

app (external e) multimetrics (string directory, string model, string reference, string out_dir) 
{
   multimetrics "-b" "1"
                "-n" "1"
                "-d" directory
                "-m" model
                "-r" reference
                "-o" out_dir;
}

string models[]   = strsplit(arg("m"), ",");
string directory  = arg("d");
string reference  = arg("r");
string out_dir    = arg("o");

foreach m in models {
   external e;
   e = multimetrics(directory, m, reference, out_dir);
}
