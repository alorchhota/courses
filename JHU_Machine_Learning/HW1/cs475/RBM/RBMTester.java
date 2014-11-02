package cs475.RBM;
import java.util.LinkedList;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

import cs475.Classify;
import cs475.CommandLineUtilities;

//NOTE: Do not modify this class
public class RBMTester {
	static public LinkedList<Option> options = new LinkedList<Option>();
	
	public static void main(String[] args) {
		String[] manditory_args = { "data"};
		createCommandLineOptions();
		CommandLineUtilities.initCommandLineParameters(args, Classify.options, manditory_args);
		
		String data_file = CommandLineUtilities.getOptionValue("data");
				
		int gd_iterations = 20;
		if (CommandLineUtilities.hasArg("gd_iterations"))
			gd_iterations = CommandLineUtilities.getOptionValueAsInt("gd_iterations");
		double gd_eta = 0.1;
		if (CommandLineUtilities.hasArg("gd_eta"))
			gd_eta = CommandLineUtilities.getOptionValueAsFloat("gd_eta");
		learning(data_file, gd_iterations, gd_eta);
	}

	public static void learning(String data_file, int gd_iterations, double gd_eta) {
		RBMParameters parameters = new RBMParameters(data_file);
		RBMEnergy rbm = new RBMEnergy(parameters, gd_iterations, gd_eta);
		rbm.learning();
		rbm.printParameters();
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
		registerOption("gd_eta", "double", true, "The constant scalar for learning rate.");
		registerOption("gd_iterations", "int", true, "The number of iterations.");
		// Other options will be added here.
	}
}


