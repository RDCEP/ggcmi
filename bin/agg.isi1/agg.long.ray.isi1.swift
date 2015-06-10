type file;

app (file o) agglongrayisi1 (string model, string gcm, string crop, string co2, string rcp, string var, string year) {
   agglongrayisi1 model gcm crop co2 rcp var year stdout = @o;
}

string models[] = strsplit(arg("models"), ",");
string gcms[]   = strsplit(arg("gcms"), ",");
string crops[]  = strsplit(arg("crops"), ",");
string co2s[]   = strsplit(arg("co2s"), ",");
string rcps[]   = strsplit(arg("rcps"), ",");
string vars[]   = strsplit(arg("vars"), ",");
string years[]  = strsplit(arg("years"), ",");

foreach m in models {
   foreach g in gcms {
      foreach c in crops {
         foreach co in co2s {
            foreach r in rcps {
               foreach v in vars {
                  foreach y in years {
                     file logfile <single_file_mapper; file = strcat("logs/", m, ".", g, ".", c, ".", co, ".", r, ".", v, ".", y, ".txt")>;
                     logfile = agglongrayisi1(m, g, c, co, r, v, y);
                  }
               }
            }
         }
      }
   }
}
