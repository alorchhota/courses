//NOTE: Do not modify this file
package cs475.RBM2;
import java.util.LinkedList;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

import cs475.Classify;
import cs475.CommandLineUtilities;
import cs475.RBM2.*;


public class RBMTester {
	static public LinkedList<Option> options = new LinkedList<Option>();
	
	public static void main(String[] args) {
		String[] manditory_args = { "data"};
		createCommandLineOptions();
		CommandLineUtilities.initCommandLineParameters(args, Classify.options, manditory_args);
		
		String data_file = CommandLineUtilities.getOptionValue("data");
				
		int num_samples = 100;
		if (CommandLineUtilities.hasArg("num_samples"))
		    num_samples = CommandLineUtilities.getOptionValueAsInt("num_samples");
		
		computeEnergy(data_file, num_samples);

	}

	public static void computeEnergy(String data_file, int num_samples) {
		//NOTE: Do not modify this function
		
		RBMParameters parameters = new RBMParameters(data_file);
		RBMEnergy rbm = new RBMEnergy(parameters, num_samples);
		for (int j=1;j<=parameters.numHiddenNodes();j++) {
			double prob = rbm.computeMarginal(j);
			System.out.println("The marginal of h_" + j + " is: " + prob);
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
		registerOption("num_samples", "int", true, "The number of samples to be generated.");
	}
}


