package cs475;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Scanner;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

public class BernoulliLikelihood {
	static public LinkedList<Option> options = new LinkedList<Option>();
	
	public static void main(String[] args) throws IOException {
		// Parse the command line.
		String[] manditory_args = { "data"};
		createCommandLineOptions();
		CommandLineUtilities.initCommandLineParameters(args, BernoulliLikelihood.options, manditory_args);
		
		String dataFile = CommandLineUtilities.getOptionValue("data");
		BernoulliLikelihood bl = new BernoulliLikelihood();
		ArrayList<Integer> data = bl.readData(dataFile);
		
		double parameter = bl.computeMaximumLikelihood(data);
		double llhood = bl.computeLogLikelihood(data, parameter);
		System.out.println("Maximum Likelihood Solution: " + Double.toString(parameter));
		System.out.println("Log-likelihood: " + Double.toString(llhood));

		
	}
	
	public static void registerOption(String option_name, String arg_name, boolean has_arg, String description) {
		OptionBuilder.withArgName(arg_name);
		OptionBuilder.hasArg(has_arg);
		OptionBuilder.withDescription(description);
		Option option = OptionBuilder.create(option_name);
		
		BernoulliLikelihood.options.add(option);		
	}
	
	private static void createCommandLineOptions() {
		registerOption("data", "String", true, "The data file to read.");
	}
	
	public double computeMaximumLikelihood(ArrayList<Integer> data) {
		// TODO: Fill in here
		// arithmatic mean is the maximum likelihood of bernoulli distribution.
		return mean(data);
		
	}
	
	public double mean(ArrayList<Integer> data){
		if(data.size() == 0)
			return 0.0;
		
		double sum = 0;
		for(Integer item : data)
			sum += item;
		return sum / data.size();
	}
	
	public double computeLogLikelihood(ArrayList<Integer> data, double parameter) {
		// log_likelihood =  sum[x * ln(u) + (1-x) * ln(1-u)]
		//				  = sum(x)*ln(u) + (N - sum(x)) * ln(1-u)
		
		int N = data.size();
		if(N==0 || parameter==0 || parameter==1)
			return 0;
		
		// 
		int s = 0;
		for(Integer item : data){
			s += item;
		}
		
		double logLikelihood = s * Math.log(parameter) + (N-s) * Math.log(1-parameter);
		return logLikelihood;
	}
	
	public ArrayList<Integer> readData(String filename) throws FileNotFoundException {
		Scanner scanner = new Scanner(new BufferedInputStream(new FileInputStream(filename)));
		ArrayList<Integer> data = new ArrayList<Integer>();
		
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			if (line.trim().length() == 0)
				   continue;
			String result = line.trim();
			int value = Integer.parseInt(result);
			data.add(value);
		}
		
		scanner.close();
		return data;
	}
}
