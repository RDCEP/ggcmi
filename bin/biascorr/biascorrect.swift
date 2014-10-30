type file;

app (file o) biascorrect(string inputfile, string reffile, string outdir) {
    biascorrect "-i" inputfile "-r" reffile "-o" outdir stdout = @o;
}

file aggfiles[] <filesys_mapper; location = "/project/joshuaelliott/ggcmi/processed/aggs/gadm0", pattern = "*">;
string reffile = "/project/joshuaelliott/ggcmi/reference/faostat/faostat.1961-2012.new.nc4";
string outdir  = "/project/joshuaelliott/ggcmi/processed/biascorr.new";

foreach f, idx in aggfiles {
    file fout <single_file_mapper; file = @strcat("./logs/bc_", idx + 1, ".out")>;
    string fn[] = @strsplit(@f, "__root__");
    fout = biascorrect(fn[1], reffile, outdir);
}
