import org.cenicana.bio.io.CollinearityParser;
import org.cenicana.bio.model.SyntenicBlock;
import java.util.List;
public class TestParser {
    public static void main(String[] args) throws Exception {
        CollinearityParser p = new CollinearityParser();
        List<SyntenicBlock> blocks = p.parse("dataset_genomico/comparativas/CC1940_vs_R570/CC1940_vs_R570.collinearity");
        int minusCount = 0;
        for (SyntenicBlock b : blocks) {
            if ("minus".equals(b.getOrientation())) minusCount++;
        }
        System.out.println("Minus blocks: " + minusCount);
        System.out.println("Total blocks: " + blocks.size());
    }
}
