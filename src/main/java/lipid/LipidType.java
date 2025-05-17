package lipid;

public enum LipidType {
    //se les asigna directamente la prioridad
    PG(1), PE(2), PI(3), PA(4), PS(5), PC(6), TG(0);

    private final int priority;

    LipidType(int priority) {
        this.priority = priority;
    }

    public int getPriority() {
        return priority;
    }
}