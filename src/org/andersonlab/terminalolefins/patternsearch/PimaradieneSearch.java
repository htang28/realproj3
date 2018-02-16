/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

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
 * of biological reactions to pull out pimaradiene-like structures.  That list of chems
 * is from the Act db from Kegg, BRENDA, and MetaCyc. The reactions have been 
 * desalted (I think).
 * 
 * See issue: Revisiting L2 terminal olefins search #1 for more explanation
 * https://github.com/UCB-BioE-Anderson-Lab/TerminalOlefins/issues/1
 * 
 * @author J. christopher Anderson
 */
public class PimaradieneSearch {
    private final String pimaradiene_smarts = "CC1(CCC2C(C1)CCC3C2(CCCC3(C)C)C)C([H])=C([H])([H])";
    private MolSearch searcher;
    
    public static void main(String[] args) throws Exception {
        //Do the search
        PimaradieneSearch search = new PimaradieneSearch();
        search.initialize();
        Set<String> hits = search.run();
        
        //Print out the hits and images
        File dir = new File("PimarSearchHits");
        dir.mkdir();
        StringBuilder sb = new StringBuilder();
        int count = 0;
        for(String inchi : hits) {
            sb.append(count).append("\t").append(inchi).append("\n");
            Molecule mol = MolImporter.importMol(inchi);
            ChemAxonUtils.savePNGImage(mol, "PimarSearchHits/" + count + ".png");
            count++;
        }
        FileUtils.writeFile(sb.toString(), "PimarSearchHits/inchis.txt");
    }

    public void initialize() throws Exception {
        //Initialize ChemAxon library (requires license file on path)
        ChemAxonUtils.license();
        
        //Construct the MolSearch query
        //From https://docs.chemaxon.com/display/jchembase/Bond+specific+search+options
        searcher = new MolSearch();
        MolSearchOptions searchOptions = new MolSearchOptions(SearchConstants.SUBSTRUCTURE);
        
        /*****    Doing Vauge Match instead of exact match *******/
        searchOptions.setVagueBondLevel(SearchConstants.VAGUE_BOND_LEVEL4);
        searcher.setSearchOptions(searchOptions);

        // queryMode = true forces string to be imported as SMARTS
        // If SMILES import needed, set queryMode = false.
        MolHandler mh1 = new MolHandler(pimaradiene_smarts, false);
        Molecule query = mh1.getMolecule();
        searcher.setQuery(query);
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
        
        //If it's not 20 carbons, ditch it
        String formula = target.getFormula();
        if(!formula.startsWith("C20")) {
            return false;
        }
        
        //If it has oxygen, ditch it
        if(formula.contains("O")) {
            return false;
        }
        
        //Do Molsearch
        searcher.setTarget(target);

        //Search all matching substructures
        //from https://www.chemaxon.com/jchem/doc/dev/java/api/chemaxon/sss/search/MolSearch.html
        int[][] hits = searcher.findAll();
        
        //Return true if the inchi contained the pattern
        if (hits == null) {
            return false;
        }
      
        return true;
    }

}
