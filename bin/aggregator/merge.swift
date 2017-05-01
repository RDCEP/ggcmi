type file;

app (file o) get_inputs (string params) {
   inputs_merge params stdout = @o;
}

app merge (string indir, string model, string crop, string clevs, string tlevs, string wlevs, string nlevs,
           string adaptation, string output) {
    merge "--indir" indir "--model" model "--crop" crop "--carbon_levels" clevs "--temp_levels" tlevs
          "--water_levels" wlevs "--nitrogen_levels" nlevs "--adaptation" adaptation "--output" output;
}

type Inputs {
    string indir;
    string model;
    string crop;
    string clevs;
    string tlevs;
    string wlevs;
    string nlevs;
    string adaptation;
    string output;
}

string params = arg("params");
file ff <"finder_merge.out">;
ff = get_inputs(params);
Inputs input[] = readData(ff);

foreach i, idx in input {
    merge(i.indir, i.model, i.crop, i.clevs, i.tlevs, i.wlevs, i.nlevs, i.adaptation, i.output);
}
