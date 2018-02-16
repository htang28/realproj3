package org.andersonlab.pathwayanalysis;

import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;
import org.andersonlab.pathwayanalysis.model.*;
import org.ucb.act.utils.ChemAxonUtils;
import org.ucb.act.utils.FileUtils;

/**
 * A Function that validates each reaction in a pathway for mass balance.
 * This is necessary to validate the stoichiometric coefficients on
 * each reaction.  This returns true if all reactions in the pathway
 * are OK.
 * 
 * @author J. Christopher Anderson
 */
public class MassValidator {
    
    public boolean run(Pathway path) throws Exception {
        //Are all the reactions balanced?
        for(Rxn rxn : path.getReactions()) {
            //Add up masses of substrates
            double subsMw = 0.0;
            for(Chem achem : rxn.getSubstrates().keySet()) {
                int stoich = rxn.getSubstrates().get(achem);
                Molecule mol = null;
                try {
                    mol = MolImporter.importMol(achem.getInchi());
                } catch(Exception err) {
                    System.out.println("Unable to parse:");
                    System.out.println(achem.getName());
                    System.out.println(achem.getInchi());
                }
                subsMw += stoich*mol.getExactMass();
            }
            
            //Add up masses of products
            double pdtsMw = 0.0;
            for (Chem achem : rxn.getProducts().keySet()) {
                int stoich = rxn.getProducts().get(achem);
                Molecule mol = null;
                try {
                    mol = MolImporter.importMol(achem.getInchi());
                } catch(Exception err) {
                    System.out.println("Unable to parse:");
                    System.out.println(achem.getName());
                    System.out.println(achem.getInchi());
                }
                pdtsMw += stoich*mol.getExactMass();
            }
            
            if(Math.abs(subsMw - pdtsMw) > 0.0001) {
                System.out.println("Balance error for\n" + rxn.toString());
                double dif = subsMw - pdtsMw;
                System.out.println("\toff by " + dif);
                return false;
            }
        }
        
        //If it survived all that, return true
        return true;
    }
    
    public static void main(String[] args) throws Exception {
        ChemAxonUtils.license();
        
        PathwayParser parser = new PathwayParser();
        String data = FileUtils.readFile("/Users/jca20n/TerminalOlefins/data/pimar.txt");
        Pathway path = parser.run(data);
        
        //Test validation
        MassValidator validator = new MassValidator();
        boolean result = validator.run(path);
        System.out.println(result);
    }
}
