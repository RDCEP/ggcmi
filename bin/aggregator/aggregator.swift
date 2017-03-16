type file;

app (file o) get_inputs () {
   inputs stdout = @o;
}

app (file o) aggregator (string indir, string crop, string lufile, string aggfile, string gsfile, string outfile) {
   aggregator "-i" indir "-c" crop "-l" lufile "-a" aggfile "-g" gsfile "-o" outfile stdout = @o;
}

type Inputs {
    string indir;
    string crop;
    string lufile;
    string agg;
    string gsfile;
    string outfile;
}

file ff <"finder.out">;
ff = get_inputs();
Inputs input[] = readData(ff);

foreach i, idx in input {
    file logfile <single_file_mapper; file = strcat("logs/log.", idx + 1, ".txt")>;
    logfile = aggregator(i.indir, i.crop, i.lufile, i.agg, i.gsfile, i.outfile);
}
