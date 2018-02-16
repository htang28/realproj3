/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.andersonlab.pathwayanalysis.model;

import org.andersonlab.pathwayanalysis.model.Chem;
import java.util.Map;

/**
 *
 * @author J. Christopher Anderson
 */
public class Rxn {
    private final Map<Chem, Integer> substrates;
    private final Map<Chem, Integer> products;

    public Rxn( Map<Chem, Integer> substrates, Map<Chem, Integer> products) {
        this.substrates = substrates;
        this.products = products;
    }


    public Map<Chem, Integer> getSubstrates() {
        return substrates;
    }

    public Map<Chem, Integer> getProducts() {
        return products;
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        Map<Chem, Integer> subs = getSubstrates();
        sb.append("Substrates: ");
        for (Chem achem : subs.keySet()) {
            int stoich = subs.get(achem);
            sb.append(stoich);
            sb.append("*");
            sb.append(achem.getName());
            sb.append(", ");
        }
        sb.append("\n");
        Map<Chem, Integer> pdts = getProducts();
        sb.append("Products: ");
        for (Chem achem : pdts.keySet()) {
            int stoich = pdts.get(achem);
            sb.append(stoich);
            sb.append("*");
            sb.append(achem.getName());
            sb.append(", ");
        }
        sb.append("\n");
        return sb.toString();
    }
}
