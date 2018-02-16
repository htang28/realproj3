package org.andersonlab.terminalolefins.patternsearch;

import chemaxon.formats.MolImporter;
import chemaxon.sss.SearchConstants;
import chemaxon.sss.search.MolSearch;
import chemaxon.sss.search.MolSearchOptions;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;
import java.io.File;
import org.andersonlab.terminalolefins.utils.ChemAxonUtils;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.andersonlab.terminalolefins.utils.FileUtils;

/**
 * This class does a substructure query against the list of products
 * of biological reactions to pull out indoxyl-like structures.  That list of chems
 * is from the Act db from Kegg, BRENDA, and MetaCyc. The reactions have been 
 * desalted (I think).
 * 
 * @author J. christopher Anderson
 */
public class IndoxylSearch {
    private final String indoxyl_smarts = "C12=C(C=CN2)C=CC=C1";
    private final String meth_indoxyl_smarts = "C12=C(C(C)=CN2)C=CC=C1";
    private MolSearch searcher;
    private MolSearch searcher2;
    
    public static void main(String[] args) throws Exception {
        //Do the search
        IndoxylSearch search = new IndoxylSearch();
        search.initialize();
        Set<String> hits = search.run();
        
        //Print out the hits and images
        File dir = new File("IndoxylHits");
        dir.mkdir();
        StringBuilder sb = new StringBuilder();
        int count = 0;
        for(String inchi : hits) {
            sb.append(count).append("\t").append(inchi).append("\n");
            Molecule mol = MolImporter.importMol(inchi);
            ChemAxonUtils.savePNGImage(mol, "IndoxylHits/" + count + ".png");
            count++;
        }
        FileUtils.writeFile(sb.toString(), "IndoxylHits/inchis.txt");
    }

    public void initialize() throws Exception {
        //Initialize ChemAxon library (requires license file on path)
        ChemAxonUtils.license();
        
        //Construct the MolSearch query
        //From https://docs.chemaxon.com/display/jchembase/Bond+specific+search+options
        searcher = new MolSearch();
        MolSearchOptions searchOptions = new MolSearchOptions(SearchConstants.SUBSTRUCTURE);
//        searchOptions.setExactBondMatching(false );
        
        /*****    Doing Vauge Match instead of exact match *******/
        searchOptions.setVagueBondLevel(SearchConstants.VAGUE_BOND_LEVEL4);
        searcher.setSearchOptions(searchOptions);

        // queryMode = true forces string to be imported as SMARTS
        // If SMILES import needed, set queryMode = false.
        MolHandler mh1 = new MolHandler(indoxyl_smarts, true);
        Molecule query = mh1.getMolecule();
        searcher.setQuery(query);
        
        
        
        //Construct the second (negative) MolSearch query
        searcher2 = new MolSearch();
        MolSearchOptions searchOptions2 = new MolSearchOptions(SearchConstants.SUBSTRUCTURE);
        searchOptions2.setVagueBondLevel(SearchConstants.VAGUE_BOND_LEVEL4);
        searcher2.setSearchOptions(searchOptions2);
        MolHandler mh2 = new MolHandler(meth_indoxyl_smarts, true);
        Molecule query2 = mh2.getMolecule();
        searcher2.setQuery(query2);
    }

    public Set<String> run() throws Exception {
        //Keep a log of failed inchis and matching inchis
        Set<String> failedInchis = new HashSet<>();
        Set<String> matchingInchis = new HashSet<>();
        
        //Parse all product molecules of reactions from Act DB
        String data = FileUtils.readFile("AllProductInchis.txt");
        String[] lines = data.split("\\r|\\r?\\n");
        for(String inchi : lines) {
            try {
                if(analyze(inchi)) {
                    matchingInchis.add(inchi);
                }
            } catch(Exception err) {
                failedInchis.add(inchi);
            }
        } 
        
        return matchingInchis;
    }
    
    private boolean analyze(String inchi) throws Exception {
        //Import the target chemical
        Molecule target  = MolImporter.importMol(inchi);
       
        //Do Molsearch for indole pattern
        searcher.setTarget(target);

        //Search all matching substructures
        //from https://www.chemaxon.com/jchem/doc/dev/java/api/chemaxon/sss/search/MolSearch.html
        int[][] hits = searcher.findAll();
        
        //Return true if the inchi contained the pattern
        if (hits == null) {
            return false;
        }
        
        //Do Second subtractive MolSearch for 3-methylindole
        searcher2.setTarget(target);

        //Search all matching substructures
        //from https://www.chemaxon.com/jchem/doc/dev/java/api/chemaxon/sss/search/MolSearch.html
        hits = searcher2.findAll();
        
        //Return true if the inchi contained the pattern
        if (hits != null) {
            return false;
        }
        
        return true;
    }

}
