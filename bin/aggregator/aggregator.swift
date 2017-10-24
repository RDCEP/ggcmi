type file;

app (file o) get_inputs (string params) {
   inputs params stdout = @o;
}

app (file o) aggregator (string indir, string crop, string lufile, string aggfile, string gsfile, string year_start, string outfile) {
   aggregator "-i" indir "-c" crop "-l" lufile "-a" aggfile "-g" gsfile "-y" year_start "-o" outfile stdout = @o;
}

type Inputs {
    string indir;
    string crop;
    string lufile;
    string agg;
    string gsfile;
    string year_start;
    string outfile;
}

string params = arg("params");
file ff <"finder.out">;
ff = get_inputs(params);
Inputs input[] = readData(ff);

foreach i, idx in input {
    file logfile <single_file_mapper; file = strcat("logs/log.", idx + 1, ".txt")>;
    logfile = aggregator(i.indir, i.crop, i.lufile, i.agg, i.gsfile, i.year_start, i.outfile);
}
