std::ofstream out("hist2d.csv");
int nx = h2->GetNbinsX();
int ny = h2->GetNbinsY();

for (int i = 1; i <= nx; ++i) {
    for (int j = 1; j <= ny; ++j) {
        double x = h2->GetXaxis()->GetBinCenter(i);
        double y = h2->GetYaxis()->GetBinCenter(j);
        double z = h2->GetBinContent(i, j);
        out << x << "," << y << "," << z << "\n";
    }
}
out.close();
