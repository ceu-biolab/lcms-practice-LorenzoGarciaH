package adduct;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MassTransformation {
    /**
     * Calculate the mass to search depending on the adduct hypothesis
     *
     * @param mz mz
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the monoisotopic mass of the experimental mass mz with the adduct @param adduct
     */
    public static Double getMonoisotopicMassFromMZ(Double mz, String adduct) {
        Double massToSearch;
        Double adductMass = getAdductMass(adduct);
        if (adductMass == null) {
            throw new IllegalArgumentException("Unknown adduct: " + adduct);
        }

        int multimer = extractMultimer(adduct);
        int charge = extractCharge(adduct);

        massToSearch = (mz + adductMass) * charge;
        return massToSearch / multimer;
    }
    private static int extractMultimer(String adduct) {
        Pattern pattern = Pattern.compile("\\[(\\d*)M.*]");
        Matcher matcher = pattern.matcher(adduct);
        if (matcher.find()) {
            String multimerStr = matcher.group(1);
            return multimerStr.isEmpty() ? 1 : Integer.parseInt(multimerStr);
        }
        return 1;
    }

    private static int extractCharge(String adduct) {
        Pattern pattern = Pattern.compile("](\\d*)([+-])");
        Matcher matcher = pattern.matcher(adduct);
        if (matcher.find()) {
            String chargeStr = matcher.group(1);
            return chargeStr.isEmpty() ? 1 : Integer.parseInt(chargeStr);
        }
        return 1;
    }
    private static Double getAdductMass(String adduct) {
        if (AdductList.MAPMZPOSITIVEADDUCTS.containsKey(adduct)) {
            return AdductList.MAPMZPOSITIVEADDUCTS.get(adduct);
        } else if (AdductList.MAPMZNEGATIVEADDUCTS.containsKey(adduct)) {
            return AdductList.MAPMZNEGATIVEADDUCTS.get(adduct);
        }
        return null;
    }

    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     *
     * @param monoisotopicMass
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return
     */
    public static Double getMZFromMonoisotopicMass(Double monoisotopicMass, String adduct) {
        Double massToSearch;
        Double adductMass = getAdductMass(adduct);
        if (adductMass == null) {
            throw new IllegalArgumentException("Unknown adduct: " + adduct);
        }

        int multimer = extractMultimer(adduct);
        int charge = extractCharge(adduct);

        massToSearch = monoisotopicMass * multimer;

        return (massToSearch - adductMass) / charge;
    }

    /**
     * Converts an m/z value to its corresponding monoisotopic mass using the adduct definition.
     *
     * @param mz     experimental mass-to-charge ratio
     * @param adduct adduct string (e.g. [M+H]+, [2M+Na]+)
     * @return monoisotopic mass of a single molecule
     */
    public static Double mzToMonoisotopicMass(Double mz, String adduct) {
        return getMonoisotopicMassFromMZ(mz, adduct);
    }

    /**
     * Converts a monoisotopic mass to the corresponding m/z value using the adduct definition.
     *
     * @param monoisotopicMass neutral monoisotopic mass of the compound
     * @param adduct           adduct string (e.g. [M+H]+, [2M+Na]+)
     * @return calculated m/z value
     */
    public static Double monoisotopicMassToMz(Double monoisotopicMass, String adduct) {
        return getMZFromMonoisotopicMass(monoisotopicMass, adduct);
    }
}
