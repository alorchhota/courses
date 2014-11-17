
// NOTE: you should not modify this code

package cs475.loopMRF;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import cs475.CommandLineUtilities;
import cs475.Classify;

public class LoopMRFTester {

	public static void main(String[] args) {
		if(args.length == 0) {
			System.out.println("provide some MRF potential data files");
			return;
		}
		String[] manditory_args = {"data"};
		createCommandLineOptions();
		CommandLineUtilities.initCommandLineParameters(args, Classify.options, manditory_args);
		String data_file = CommandLineUtilities.getOptionValue("data");
		int iterations = 50;
		if (CommandLineUtilities.hasArg("iterations"))
			iterations = CommandLineUtilities.getOptionValueAsInt("iterations");
		findPotentials(data_file, iterations);
	}

	public static void findPotentials(String data_file, int iterations) {
		LoopMRFPotentials p = new LoopMRFPotentials(data_file);
		LoopyBP lbp = new LoopyBP(p, iterations);
		for(int i=1; i<=p.loopLength(); i++) {
			double[] marginal = lbp.marginalProbability(i);
			if(marginal.length-1 != p.numXValues())		// take off 1 for 0 index which is not used
				throw new RuntimeException("length of probability distribution is incorrect: " + marginal.length);
			System.out.println("marginal probability distribution for node " + i + " is:");
			double sum = 0.0;
			for(int k=1; k<=p.numXValues(); k++) {
				if(marginal[k] < 0.0 || marginal[k] > 1.0)
					throw new RuntimeException("illegal probability for x_" + i);
				System.out.println("\tP(x = " + k + ") = " + marginal[k]);
				sum += marginal[k];
			}
			double err = 1e-5;
			if(sum < 1.0-err || sum > 1.0+err)
				throw new RuntimeException("marginal probability distribution for x_" + i + " doesn't sum to 1");
		}
	}

	public static void registerOption(String option_name, String arg_name, boolean has_arg, String description) {
		OptionBuilder.withArgName(arg_name);
		OptionBuilder.hasArg(has_arg);
		OptionBuilder.withDescription(description);
		Option option = OptionBuilder.create(option_name);
		
		Classify.options.add(option);		
	}
	
	private static void createCommandLineOptions() {
		registerOption("data", "String", true, "The data to use.");
		registerOption("iterations", "int", true, "The number of iterations. default: 50");
		// Other options will be added here.
	}
}

