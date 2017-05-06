type file;

app (file o) get_inputs (string params) {
   inputs params stdout = @o;
}

app aggregator (string indir, string crop, string lufile, string aggfile, string gsfile, string co2,
                string temperature, string precip, string nitrogen, string adaptation, string outfile) {
   aggregator "--indir" indir "--crop" crop "--lufile" lufile "--agg" aggfile "--gsfile" gsfile "--co2" co2
              "--temperature" temperature "--precipitation" precip "--nitrogen" nitrogen "--adaptation" adaptation
              "--outfiles" outfile;
}

type Inputs {
    string indir;
    string crop;
    string lufile;
    string agg;
    string gsfile;
    string co2;
    string temperature;
    string precip;
    string nitrogen;
    string adaptation;
    string outfile;
}

string params = arg("params");
file ff <"finder.out">;
ff = get_inputs(params);
Inputs input[] = readData(ff);

foreach i, idx in input {
    aggregator(i.indir, i.crop, i.lufile, i.agg, i.gsfile, i.co2, i.temperature, i.precip, i.nitrogen, i.adaptation, i.outfile);
}
