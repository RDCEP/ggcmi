type file;

app (file o) cleanlongisi1 (string model, string gcm, string crop, string co2, string rcp, string var) {
   cleanlongisi1 model gcm crop co2 rcp var stdout = @o;
}

string models[] = strsplit(arg("models"), ",");
string gcms[]   = strsplit(arg("gcms"), ",");
string crops[]  = strsplit(arg("crops"), ",");
string co2s[]   = strsplit(arg("co2s"), ",");
string rcps[]   = strsplit(arg("rcps"), ",");
string vars[]   = strsplit(arg("vars"), ",");

foreach m in models {
   foreach g in gcms {
      foreach c in crops {
         foreach co in co2s {
            foreach r in rcps {
               foreach v in vars {
                  file logfile <single_file_mapper; file = strcat("logs/", m, ".", g, ".", c, ".", co, ".", r, ".", v, ".txt")>;
                  logfile = cleanlongisi1(m, g, c, co, r, v);
               }
            }
         }
      }
   }
}
