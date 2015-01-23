type file;

app (file o) plotisi1 (string plot, string crop, string co2) {
   plotisi1 plot crop co2 stdout = @o;
}

string plots[] = strsplit(arg("plots"), ",");
string crops[] = strsplit(arg("crops"), ",");
string co2s[]  = strsplit(arg("co2s"),  ",");

foreach p in plots {
   foreach c in crops {
      foreach co in co2s {
         file logfile <single_file_mapper; file=strcat("logs/", p, ".", c, ".", co, ".txt")>;
         logfile = plotisi1(p, c, co);
      }
   }
}
