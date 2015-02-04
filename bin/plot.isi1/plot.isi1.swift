type file;

app (file o) plotisi1 (string plot, string var, string crop) {
   plotisi1 plot var crop stdout = @o;
}

string plots[] = strsplit(arg("plots"), ",");
string vars[]  = strsplit(arg("vars"),  ",");
string crops[] = strsplit(arg("crops"), ",");

foreach p in plots {
   foreach v in vars {
      foreach c in crops {
         file logfile <single_file_mapper; file=strcat("logs/", p, ".", v, ".", c, ".txt")>;
         logfile = plotisi1(p, v, c);
      }
   }
}
