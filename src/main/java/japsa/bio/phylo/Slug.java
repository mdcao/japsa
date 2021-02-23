package japsa.bio.phylo;

import java.text.Normalizer;
import java.text.Normalizer.Form;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Slug{
	private static final Pattern NONLATIN = Pattern.compile("[^\\w-]");
	private static final Pattern WHITESPACE = Pattern.compile("[\\s_]");
	private static String toSlug(String input) {
		return toSlug(input, "_");
	}
	private static String toSlug(String input, int k) {
		return toSlug(input, k, "_");
	}
	public static String toSlug(String input, String replace) {
	 input = input.replaceAll(", complete sequence", "").replaceAll("-","");
	 
	  Matcher matcher = WHITESPACE.matcher(input);
	 String nowhitespace = WHITESPACE.matcher(input).replaceAll(replace);
	  String normalized = Normalizer.normalize(nowhitespace, Form.NFD);
	  String slug = NONLATIN.matcher(normalized).replaceAll("");
	  return slug.toLowerCase(Locale.ENGLISH);
	}
	public static String toSlug(String input, int k, String replace) {
		  Matcher matcher = WHITESPACE.matcher(input);
		  boolean find = matcher.find();
		  for(int i=1;  i < k; i++) {
			  find = matcher.find();
		  }
		  if(find){
			  input = input.substring(0,matcher.start());
		  }
		 String nowhitespace = WHITESPACE.matcher(input).replaceAll(replace);
		  String normalized = Normalizer.normalize(nowhitespace, Form.NFD);
		  String slug = NONLATIN.matcher(normalized).replaceAll("");
		  return slug.toLowerCase(Locale.ENGLISH);
		}
	
}