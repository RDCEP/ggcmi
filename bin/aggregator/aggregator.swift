type file;

app (file o) aggregator (string directory, string model, string weather, string crop, string landuse_ir, 
                             string landuse_rf, string aggregation_mask, string timestamp, string growing_dir,
                             string output_dir) 
{
   aggregator "-b" "1" 
              "-n" "1"
              "-d" directory
              "-m" model
              "-w" weather
              "-c" crop
              "-i" landuse_ir
              "-r" landuse_rf
              "-a" aggregation_mask
              "-t" timestamp
              "-g" growing_dir
              "-o" output_dir
              stdout=@o;
}

string directory  = arg("d");
string models[]   = strsplit(arg("m"), ",");
string weather[]  = strsplit(arg("w"), ",");
string crop[]     = strsplit(arg("c"), ",");
string landuse_ir = arg("i");
string landuse_rf = arg("r");
string agg_mask[] = strsplit(arg("a"), ",");
string timestamp  = arg("t");
string grow_dir   = arg("g");
string out_dir[]  = strsplit(arg("o"), ",");

foreach m in models {
   foreach w in weather {
      foreach c in crop {
         foreach a,idx in agg_mask {
            file logfile <single_file_mapper; file=strcat("logs/", m, ".", w, ".", c, ".mask_", idx, ".txt")>;
            logfile = aggregator(directory, m, w, c, landuse_ir, landuse_rf, agg_mask[idx], timestamp, grow_dir, out_dir[idx]);
         }
      }
   }
}
