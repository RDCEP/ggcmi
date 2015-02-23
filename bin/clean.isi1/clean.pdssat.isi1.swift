type file;

app (file o) cleanpdssatisi1 (string gcm, string crop, string co2, string rcp) {
   cleanpdssatisi1 gcm crop co2 rcp stdout = @o;
}

string gcms[]  = strsplit(arg("gcms"), ",");
string crops[] = strsplit(arg("crops"), ",");
string co2s[]  = strsplit(arg("co2s"), ",");
string rcps[]  = strsplit(arg("rcps"), ",");

foreach g in gcms {
   foreach c in crops {
      foreach co in co2s {
         foreach r in rcps {
            file logfile <single_file_mapper; file = strcat("logs/", g, ".", c, ".", co, ".", r, ".txt")>;
            logfile = cleanpdssatisi1(g, c, co, r);
         }
      }
   }
}
