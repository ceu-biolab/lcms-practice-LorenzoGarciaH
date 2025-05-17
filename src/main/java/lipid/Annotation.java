package lipid;

import adduct.Adduct;
import adduct.AdductList;
import adduct.MassTransformation;

import java.util.*;

/**
 * Class to represent the annotation over a lipid
 */
public class Annotation {

    private final Lipid lipid;
    private final double mz;
    private final double intensity; // intensity of the most abundant peak in the groupedPeaks
    private final double rtMin;
    private final IoniationMode ionizationMode;
    private String adduct; // !!TODO The adduct will be detected based on the groupedSignals
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;


    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode) {
        this(lipid, mz, intensity, retentionTime, ionizationMode, Collections.emptySet());
    }

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     * @param groupedSignals
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode, Set<Peak> groupedSignals) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.ionizationMode = ionizationMode;
        // !!TODO This set should be sorted according to help the program to deisotope the signals plus detect the adduct
        this.groupedSignals = new TreeSet<>(Comparator.comparingDouble(Peak::getMz));
        this.groupedSignals.addAll(groupedSignals); // add all the signals if not is empty
        this.score = 0;
        this.totalScoresApplied = 0;
        this.adduct = detectAdduct();
    }

    public Lipid getLipid() {
        return lipid;
    }

    public double getMz() {
        return mz;
    }

    public double getRtMin() {
        return rtMin;
    }

    public String getAdduct() {
        return adduct;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    public double getIntensity() {
        return intensity;
    }

    public IoniationMode getIonizationMode() {
        return ionizationMode;
    }

    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }


    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    // !CHECK Take into account that the score should be normalized between -1 and 1
    public void addScore(int delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    /**
     * @return The normalized score between 0 and 1 that consists on the final number divided into the times that the rule
     * has been applied.
     */
    public double getNormalizedScore() {
        return (double) this.score / this.totalScoresApplied;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation)) return false;
        Annotation that = (Annotation) o;
        return Double.compare(that.mz, mz) == 0 &&
                Double.compare(that.rtMin, rtMin) == 0 &&
                Objects.equals(lipid, that.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format("Annotation(%s, mz=%.4f, RT=%.2f, adduct=%s, intensity=%.1f, score=%d)",
                lipid.getName(), mz, rtMin, adduct, intensity, score);
    }

    public String detectAdduct() {
        double mzTolerance = 0.1;
        if (groupedSignals == null || groupedSignals.size() < 2) {
            System.out.println("detectAdduct: Not enough signals (" +
                    (groupedSignals == null ? 0 : groupedSignals.size()) + ")");
            return null;
        }

        Map<String, Double> adductMap = (getIonizationMode() == IoniationMode.POSITIVE)
                ? AdductList.MAPMZPOSITIVEADDUCTS
                : AdductList.MAPMZNEGATIVEADDUCTS;

        double observedMz = this.getMz();
        System.out.println("detectAdduct: observedMz = " + observedMz + ", mode = " + getIonizationMode());

        // for each adduct that could detect observedMz:
        for (String candidateAdduct : adductMap.keySet()) {
            System.out.println("  Probando candidateAdduct: " + candidateAdduct);
            try {
                // monoisotopic mass per adducto for observedMz
                double monoisotopicMass = MassTransformation.getMonoisotopicMassFromMZ(observedMz, candidateAdduct);
                System.out.println("    monoisotopicMass = " + monoisotopicMass);

                for (Peak otherPeak : groupedSignals) {
                    System.out.println("    Comparando con otherPeak: " + otherPeak);
                    if (Math.abs(otherPeak.getMz() - observedMz) <= mzTolerance) {
                        // Same peak so skip
                        System.out.println("      Skip: same as observedMz");
                        continue;
                    }

                    for (String secondAdduct : adductMap.keySet()) {
                        double expectedMz = MassTransformation
                                .getMZFromMonoisotopicMass(monoisotopicMass, secondAdduct);
                        double diff = Math.abs(expectedMz - otherPeak.getMz());
                        System.out.println("      secondAdduct=" + secondAdduct + ", expectedMz=" + expectedMz + ", observed=" + otherPeak.getMz() + ", diff=" + diff);
                        if (diff <= mzTolerance) {
                            System.out.println("    DETECTED adduct: " + candidateAdduct + " (via " + secondAdduct + ")");
                            return candidateAdduct;
                        }
                    }
                }

            } catch (IllegalArgumentException e) {

                System.out.println(" redundant candidate Adduct : " + candidateAdduct);
            }
        }

        System.out.println("detectAdduct: no adducts are accepted");
        return null;

    }

    // !!TODO Detect the adduct with an algorithm or with drools, up to the  user.
}
