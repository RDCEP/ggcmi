type file;

app (file o) metricsisi1 (string model, string gcm, string crop, string co2) {
   metricsisi1 model gcm crop co2 stdout = @o;
}

string models[] = strsplit(arg("models"), ",");
string gcms[]   = strsplit(arg("gcms"), ",");
string crops[]  = strsplit(arg("crops"), ",");
string co2s[]   = strsplit(arg("co2s"), ",");

foreach m in models {
   foreach g in gcms {
      foreach c in crops {
         foreach co in co2s {
            file logfile <single_file_mapper; file=strcat("logs/", m, ".", g, ".", c, ".", co, ".txt")>;
            logfile = metricsisi1(m, g, c, co);
         }
      }
   }
}
