type file;

app (file o) plotisi1 (string plot, string crop) {
   plotisi1 plot crop stdout = @o;
}

string plots[] = strsplit(arg("plots"), ",");
string crops[] = strsplit(arg("crops"), ",");

foreach p in plots {
   foreach c in crops {
      file logfile <single_file_mapper; file=strcat("logs/", p, ".", c, ".txt")>;
      logfile = plotisi1(p, c);
   }
}
