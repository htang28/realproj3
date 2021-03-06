package org.ucb.pathwayanalysis;

import org.ucb.pathwayanalysis.model.Chem;
import org.ucb.pathwayanalysis.model.Rxn;
import org.ucb.pathwayanalysis.model.Pathway;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.ucb.act.utils.ChemAxonUtils;
import org.ucb.act.utils.FileUtils;

/**
 * A Function that inputs a Pathway, performs flux balance analysis, and outputs
 * the flux vector. You may need to add the names of additional cofactors
 * included in your pathway.
 *
 * @author J. Christopher Anderson
 */
public class PathwayBalancer {

    private final Set<String> cofactors;

    public PathwayBalancer() {
        //A universal list of standard-named cofactors, you may need to add more
        cofactors = new HashSet<>();
        cofactors.add("ATP");
        cofactors.add("ADP");
        cofactors.add("H+");
        cofactors.add("H2O");
        cofactors.add("NAD+");
        cofactors.add("NADH");
        cofactors.add("NADP+");
        cofactors.add("NADPH");
        cofactors.add("Pi");
        cofactors.add("FAD");
        cofactors.add("FADH2");
        cofactors.add("CO2");
        cofactors.add("CoA");
    }

    public Map<Chem, Integer> run(Pathway path) throws Exception {
        //Pull out all non-cofactor, non input chems and fix them to an index
        Map<Chem, Integer> intermedToIndex = new HashMap<>();
        int chemCount = 0;

        for (Rxn rxn : path.getReactions()) {
            Map<Chem, Integer> allChems = rxn.getSubstrates();
            Map<Chem, Integer> pdtToStoichiometry = rxn.getProducts();

            //Add all the products to allChems (which already has the substrates)
            for (Chem achem : pdtToStoichiometry.keySet()) {
                allChems.put(achem, 0);
            }

            //For all chems (pdts or subs)
            for (Chem achem : allChems.keySet()) {
                //Ignore the cofactors
                if (cofactors.contains(achem.getName())) {
                    continue;
                }

                //Ignore the inputs
                if (path.getInputs().contains(achem)) {
                    continue;
                }

                //If it's already in the list, ignore it
                if (!intermedToIndex.containsKey(achem)) {
                    intermedToIndex.put(achem, chemCount);
                    chemCount++;
                }
            }
        }

        //See which is bigger, chemCount or rxnCount
        int rxnCount = path.getReactions().size();
        int n = Math.max(rxnCount, chemCount);

        //Construct the matrix and objective function, and zero all out
        double[][] stoichMatrix = new double[n][n];
        double[] objectiveFunc = new double[n];
        for (int i = 0; i < n; i++) {
            objectiveFunc[i] = 0;
            for (int j = 0; j < n; j++) {
                stoichMatrix[i][j] = 0;
            }
        }

        //Populate the indices of the objective function that are outputs as 1
        for (Chem outputChem : path.getOutputs()) {
            int index = intermedToIndex.get(outputChem);
            objectiveFunc[index] = 1;
        }

        //Populate the stoich matrix with all rxn stoich indices
        for (int j = 0; j < path.getReactions().size(); j++) {
            Rxn rxn = path.getReactions().get(j);

            //Put in each substrate
            for (Chem achem : rxn.getSubstrates().keySet()) {
                int stoich = rxn.getSubstrates().get(achem);
                Integer index = intermedToIndex.get(achem);
                if (index == null) {
                    continue;
                }
                stoichMatrix[index][j] = -1 * stoich;
            }

            //Put in each product
            for (Chem achem : rxn.getProducts().keySet()) {
                int stoich = rxn.getProducts().get(achem);
                Integer index = intermedToIndex.get(achem);
                if (index == null) {
                    continue;
                }
                stoichMatrix[index][j] = stoich;
            }
        }

        //Solve to get coefficients on each reaction
        SolveLinearEquation sle = new SolveLinearEquation();
        double[] result = sle.run(stoichMatrix, objectiveFunc);

        //Find the multiplier that converts the coefficients to integers
        int multiplier = 1;
        while (true) {
            double total = 0;
            for (double coeff : result) {
                total += coeff * multiplier;
            }
            double offset = total - Math.floor(total);
            if (offset < 0.001) {
                break;
            }
            multiplier++;
            if (multiplier > 20) {
                System.err.println("Unable to resolve solution to integers");
                throw new Exception();
            }
        }


        //Reinterpret the coefficients as integers
        int[] rxnCoeffs = new int[result.length];
        for (int i = 0; i < rxnCoeffs.length; i++) {
            rxnCoeffs[i] = 0;
        }
        for (int i = 0; i < result.length; i++) {
            double coeff = result[i];
            int newcoeff = (int) Math.round(multiplier * coeff);
            rxnCoeffs[i] = newcoeff;
        }
        
        //Construct the Map that will be returned by this method
        Map<Chem, Integer> bal = new HashMap<>();

        //Calculate the new balance
        for (int i = 0; i < path.getReactions().size(); i++) {
            Rxn rxn = path.getReactions().get(i);
            int rxnCoeff = rxnCoeffs[i];

            //Handle the substrates
            for (Chem achem : rxn.getSubstrates().keySet()) {
                int stoich = rxn.getSubstrates().get(achem);
                int flux = -1 * stoich * rxnCoeff;
                Integer existing = bal.get(achem);
                if (existing == null) {
                    existing = 0;
                }
                existing += flux;
                bal.put(achem, existing);
            }

            //Handle the products
            for (Chem achem : rxn.getProducts().keySet()) {
                int stoich = rxn.getProducts().get(achem);
                int flux = stoich * rxnCoeff;
                Integer existing = bal.get(achem);
                if (existing == null) {
                    existing = 0;
                }
                existing += flux;
                bal.put(achem, existing);
            }
        }

        //Delete all the zero balances
        Set<Chem> tossers = new HashSet<>();
        for (Chem achem : bal.keySet()) {
            int val = bal.get(achem);
            if (val == 0) {
                tossers.add(achem);
            }
        }

        for (Chem toss : tossers) {
            bal.remove(toss);
        }

        return bal;
    }

    public static void main(String[] args) throws Exception {
        //Parse the pathway
        PathwayParser parser = new PathwayParser();

        //Note:  You need to replace the path to point to the file on your system
        String data = FileUtils.readFile("/Users/jca20n/TerminalOlefins/data/ethanol.txt");
//        String data = FileUtils.readFile("/Users/jcaucb/Documents/TerminalOlefins/data/chorismate.txt");
//        String data = FileUtils.readFile("/Users/jcaucb/Documents/TerminalOlefins/data/pimar.txt");
//        String data = FileUtils.readFile("/Users/jcaucb/Documents/TerminalOlefins/data/mevalonate.txt");
//        String data = FileUtils.readFile("/Users/jcaucb/Documents/TerminalOlefins/data/galact.txt");
//        String data = FileUtils.readFile("/Users/jca20n/TerminalOlefins/data/abstract_example.txt");
        Pathway path = parser.run(data);
        System.out.println("Examining pathway:\n\n" + path.getName());

        //Validate the pathway
        ChemAxonUtils.license();
        MassValidator validator = new MassValidator();
        boolean result = validator.run(path);
        System.out.println("\nDoes the pathway pass mass balance validation on each reaction?");
        System.out.println(result);

        //Balance the pathway
        PathwayBalancer balancer = new PathwayBalancer();
        Map<Chem, Integer> bal = balancer.run(path);
        System.out.println("\nBalance is:\n");
        for (Chem achem : bal.keySet()) {
            System.out.println(achem.getName() + " : " + bal.get(achem));
        }
    }
}
