// V6 harness: SciPhy's own tree likelihood (azwaans/SciPhy, SciPhyTreeLikelihood, unmodified source
// compiled against BEAST 2.7.7) evaluated on trees, parameters and reads written by
// 06_validate_sciphy.py. One SciPhyTreeLikelihood per tape (SciPhy's convention: one alignment = one
// tape, columns = sites). Strict clock at the tape's rate, origin = 1, tips at height 0.
//
// Input file:  line 1  "n k N M"
//              line 2  newick with taxa t0..t{n-1}, branch lengths in experiment units
//              line 3  M insert probabilities (space separated; insert i uses entry i-1)
//              then per tape: one line with the clock rate, then n lines "s1,...,sN" (0 = unedited)
// Output:      one line per tape, the log-likelihood (nats) to 17 significant digits.
//
// usage: java -cp BEAST.base.jar:<sciphy classes>:. V6Harness in1.txt [in2.txt ...]

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import sciphy.evolution.likelihood.SciPhyTreeLikelihood;
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class V6Harness {

    public static void main(String[] args) throws Exception {
        // outside BEAST's launcher the service registry is empty: register BEAST.base's data types etc.
        beast.pkgmgmt.BEASTClassLoader.addServices(System.getProperty("beast.base.versionxml"));
        for (String path : args) {
            run(path);
        }
    }

    static void run(String path) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(path));
        String[] h = br.readLine().trim().split("\\s+");
        int n = Integer.parseInt(h[0]), k = Integer.parseInt(h[1]), N = Integer.parseInt(h[2]), M = Integer.parseInt(h[3]);
        String newick = br.readLine().trim();
        String xi = br.readLine().trim();

        StringBuilder fr = new StringBuilder("1.0");
        for (int i = 0; i < M; i++) fr.append(" 0.0");

        PrintWriter out = new PrintWriter(new java.io.FileWriter(path + ".sciphy_ll"));
        for (int z = 0; z < k; z++) {
            String rate = br.readLine().trim();
            List<Sequence> seqs = new ArrayList<>();
            for (int c = 0; c < n; c++) {
                seqs.add(new Sequence("t" + c, br.readLine().trim()));
            }
            Alignment aln = new Alignment(seqs, "integer");
            Tree tree = new TreeParser(aln, newick);

            RealParameter ep = new RealParameter();
            ep.initByName("value", xi, "dimension", M);
            RealParameter fp = new RealParameter();
            fp.initByName("value", fr.toString(), "dimension", M + 1);
            Frequencies freqs = new Frequencies();
            freqs.initByName("frequencies", fp, "estimate", false);
            SciPhySubstitutionModel sm = new SciPhySubstitutionModel();
            sm.initByName("editProbabilities", ep, "frequencies", freqs);
            SiteModel site = new SiteModel();
            site.initByName("substModel", sm, "mutationRate", new RealParameter("1.0"));
            StrictClockModel clock = new StrictClockModel();
            clock.initByName("clock.rate", new RealParameter(rate));

            SciPhyTreeLikelihood L = new SciPhyTreeLikelihood();
            L.initByName("data", aln, "tree", tree, "siteModel", site, "branchRateModel", clock,
                    "origin", new RealParameter("1.0"), "arrayLength", new IntegerParameter(Integer.toString(N)),
                    "useScaling", true);
            out.printf("%.17g%n", L.calculateLogP());
        }
        out.close();
        br.close();
    }
}
