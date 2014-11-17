type file;

app (file o) get_inputs () {
   inputs stdout = @o;  
}

app rescaler(string irfile, string rffile, string bcfile, string mkfile, string wtfile, string crmthd, string agglvl, string vname, string outfile) {
   rescaler "-i" irfile "-r" rffile "-b" bcfile "-m" mkfile "-w" wtfile "-c" crmthd "-a" agglvl "-v" vname "-o" outfile;
}

type Inputs {
   string irfile;
   string rffile;
   string bcfile;
   string wtfile;
   string vname;
   string outfile;
}

string mkfile = arg("mkfile");
string crmthd = arg("crmthd");
string agglvl = arg("agglvl");

file ff <"finder.out">;
ff = get_inputs();
Inputs inp[] = readData(ff);

foreach i in inp {
   rescaler(i.irfile, i.rffile, i.bcfile, mkfile, i.wtfile, crmthd, agglvl, i.vname, i.outfile);
}
