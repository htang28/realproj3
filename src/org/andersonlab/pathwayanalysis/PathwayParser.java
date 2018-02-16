
package org.andersonlab.pathwayanalysis;

import java.util.ArrayList;
import org.andersonlab.pathwayanalysis.model.Rxn;
import org.andersonlab.pathwayanalysis.model.Pathway;
import org.andersonlab.pathwayanalysis.model.Chem;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.ucb.act.utils.FileUtils;

/**
 *
 * @author J. Christopher Anderson
 */
public class PathwayParser {
    
    public Pathway run(String data) throws Exception {
        data = data.replaceAll("\"", "");
        String[] regions = data.split("@");
        Map<String,String> labelToText = new HashMap<>();
        for(int i=0; i<regions.length; i++) {
            if(regions[i].isEmpty()) {
                continue;
            }
            
            int firstColon = regions[i].indexOf(":");
            String label = regions[i].substring(0, firstColon);
            labelToText.put(label, regions[i].substring(firstColon + 1).trim());
        }
        
        //Handle the chemicals
        Map<String, Chem> chemicals = new HashMap<>();
        try {
            String chemData = labelToText.get("chemicals");
            String[] lines = chemData.split("\\r|\\r?\\n");
            for(String line : lines) {
                if(line.isEmpty()) {
                    continue;
                }
                String[] tabs = line.split("\t");

                if(tabs.length != 2) {
                    System.err.println("Error parsing text on:\n" + line);
                    throw new Exception();
                }

                String name = tabs[0];
                String inchi = tabs[1];
                Chem chem = new Chem(name, inchi);
                chemicals.put(name, chem);
            }
        } catch(Exception err) {
            System.err.println("Error processing the chemicals");
            throw err;
        }
        
        //Handle the inputs
        Set<Chem> inputs = new HashSet<>();
        try {
            String inputsData = labelToText.get("inputs");
            String[] lines = inputsData.split("\\r|\\r?\\n");
            for(String line : lines) {
                Chem achem = chemicals.get(line);
                if(achem == null) {
                    System.err.println("Error findig input chem: " + line);
                    throw new Exception();
                }
                inputs.add(achem);
            }
        } catch(Exception err) {
            System.err.println("Error processing the inputs");
            throw err;
        }
        
        //Handle the outputs
        Set<Chem> outputs = new HashSet<>();
        try {
            String outputsData = labelToText.get("outputs");
            String[] lines = outputsData.split("\\r|\\r?\\n");
            for(String line : lines) {
                Chem achem = chemicals.get(line);
                if(achem == null) {
                    System.err.println("Error findig output chem: " + line);
                    throw new Exception();
                }
                outputs.add(achem);
            }
        } catch(Exception err) {
            System.err.println("Error processing the outputs");
            throw err;
        }
        
        //Handle the name
        String name = labelToText.get("name");
        if(name == null) {
            System.err.println("Could not parse name");
            throw new Exception();
        }
        
        //Handle the reactions
        List<Rxn> reactions = new ArrayList<>();
        try {
            String rxnData = labelToText.get("reactions");
            String[] lines = rxnData.split("\\r|\\r?\\n");
            for(String line : lines) {
                if(line.isEmpty()) {
                    continue;
                }
                if(!line.contains(" --> ")) {
                    continue;
                }
                
                //Split the substrates and products
                String[] sides = line.split("\\s+-->\\s+");
                String subsString = sides[0];
                String pdtsString = sides[1];
                
                //Pull out the stoichiometry and chemical reference
                Map<Chem, Integer> substrates = extractChems(subsString, chemicals);
                Map<Chem, Integer> products = extractChems(pdtsString, chemicals);
                
                //Construct the reaction
                Rxn rxn = new Rxn(substrates, products);
                reactions.add(rxn);
            }
        } catch(Exception err) {
            System.err.println("Error processing the reactions");
            throw err;
        }
        
        Pathway pathway = new Pathway(name, inputs, outputs, reactions, chemicals);
        return pathway;
    }
    
    private Map<Chem, Integer> extractChems(String chemString, Map<String, Chem> chemicals) throws Exception {
        Map<Chem, Integer> out = new HashMap<>();
        String[] chemTokens = chemString.split("\\s+[+]\\s+");
        for(String token : chemTokens) {
            int stoich = 1;
            String name = null;
            
            int numterms = token.split("[0-9]+\\s").length;
            if( numterms == 2) {
                String[] stoichAndName = token.split("\\s");
                stoich = Integer.parseInt(stoichAndName[0].trim());
                name = stoichAndName[1];
            } else if( numterms == 1) {
                name = token.trim();
            } else {
                System.err.println("Error parsing chemical token: " + token);
                throw new Exception();
            }
         
            Chem achem = chemicals.get(name);
            if (achem == null) {
                System.err.println("Error extracting rxn chem: " + name);
                throw new Exception();
            }
            out.put(achem, stoich);
        }
        return out;
    }

    public static void main(String[] args) throws Exception {
        PathwayParser parser = new PathwayParser();
        String data = FileUtils.readFile("/Users/jca20n/TerminalOlefins/data/glycolysis.txt");
//        System.out.println(data);
        Pathway path = parser.run(data);
        
        System.out.println(path.getName());
        System.out.println(path.getReactions().size() == 10);
    }


}
