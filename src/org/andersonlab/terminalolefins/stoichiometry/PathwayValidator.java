/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.andersonlab.terminalolefins.stoichiometry;

import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.andersonlab.terminalolefins.stoichiometry.model.Chem;
import org.andersonlab.terminalolefins.stoichiometry.model.Pathway;
import org.andersonlab.terminalolefins.stoichiometry.model.Rxn;
import org.andersonlab.terminalolefins.utils.ChemAxonUtils;
import org.andersonlab.terminalolefins.utils.FileUtils;

/**
 *
 * @author J. Christopher Anderson
 */
public class PathwayValidator {
    
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
        PathwayValidator validator = new PathwayValidator();
        boolean result = validator.run(path);
        System.out.println(result);
    }
}
