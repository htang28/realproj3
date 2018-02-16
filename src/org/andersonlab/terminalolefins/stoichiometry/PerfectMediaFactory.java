package org.andersonlab.terminalolefins.stoichiometry;

import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.andersonlab.terminalolefins.utils.FileUtils;

/**
 *
 * @author J. Christopher Anderson
 */
public class PerfectMediaFactory {
    
    List<Chem> mediaChems;
    Map<Integer, Double> atnoToComposition;
    
    class Chem {
        String inchi;
        String name;
        String desc;
        double pricePerTon;
    }
    
    public void run() throws Exception {
        parseChemicals();
        parseTargetComposition();
        
        //Lock indicies for each atom
        int[] atomicNumbers = new int[atnoToComposition.size()];
        int x = 0;
        for(Integer atno : atnoToComposition.keySet()) {
            atomicNumbers[x] = atno;
            x++;
        }
        
        //Construct stoichiometry matrix
        double[][] stoichMatrix = new double[atomicNumbers.length][mediaChems.size()];
        for(int i=0; i<stoichMatrix.length; i++) {
            for(int j=0; j<stoichMatrix[i].length; j++) {
                stoichMatrix[i][j] = 0.0;
            }
        }
        
        for(int mc=0; mc<mediaChems.size(); mc++) {
            Chem achem = mediaChems.get(mc);
            Molecule amol = MolImporter.importMol(achem.inchi);
            for(int at=0; at<atomicNumbers.length; at++) {
                int count = amol.getAtomCount(atomicNumbers[at]);
                stoichMatrix[at][mc] = count;
            }
        }
        
        //Construct the target composition
        double[] composition = new double[atomicNumbers.length];
        for(int i=0; i<atomicNumbers.length; i++) {
            composition[i] = atnoToComposition.get(atomicNumbers[i]);
        }
        
        //Fill in a square version of the matrix with zeros
        double[][] squareMatrix = new double[composition.length][composition.length];
        for(int i=0; i<squareMatrix.length; i++) {
            for(int j=0; j<squareMatrix[i].length; j++) {
                squareMatrix[i][j] = 0.0;
            }
        }
        for(int i=0; i<stoichMatrix.length; i++) {
            for(int j=0; j<stoichMatrix[i].length; j++) {
                squareMatrix[i][j] = stoichMatrix[i][j];
            }
        }
        
        
        //Solve the stoichiometry matrix to match the composition
        double[] soln = Solve_Linear_Equation.solve(squareMatrix, composition);
        
        
        
        System.out.println("done");
    }
    
    private void parseChemicals() throws Exception {
        mediaChems = new ArrayList<>();
        String data = FileUtils.readFile("/Users/jca20n/TerminalOlefins/data/media_chems.txt");
        String[] lines = data.split("\\r|\\r?\\n");
        for(int i=1; i<lines.length; i++) {
            String line = lines[i];
            String[] tabs = line.split("\t");
            Chem achem = new Chem();
            achem.name = tabs[0];
            achem.inchi = tabs[1];
            achem.desc = tabs[2];
            achem.pricePerTon = Double.parseDouble(tabs[3]);
            mediaChems.add(achem);
        }
    }

    /**
     * Data from dry mass weights of E. coli for common elements:
     * http://www.gatewaycoalition.org/files/hidden/react/ch2/2_1f.htm	
     * 
     * H and O values are excluded because they will already
     * be in excess in water and thus are unconstrained
     * 
     * @throws Exception 
     */
    private void parseTargetComposition() throws Exception {
        atnoToComposition = new HashMap();
        String data = FileUtils.readFile("/Users/jca20n/TerminalOlefins/data/dry_mass_weights.txt");
        String[] lines = data.split("\\r|\\r?\\n");
        for(int i=2; i<lines.length; i++) {
            String line = lines[i];
            String[] tabs = line.split("\t");
            int atno = Integer.parseInt(tabs[1]);
            Double composition = Double.parseDouble(tabs[2]);
            atnoToComposition.put(atno, composition);
        }
        System.out.println("There are "  + atnoToComposition.size() +  " atoms");
    }

    public static void main(String[] args) throws Exception {
        PerfectMediaFactory pfm = new PerfectMediaFactory();
        pfm.run();
        System.out.println("done");
    }
}
