# Copyright 2017 Anurag Priyam, Queen Mary University of London
#
# Helper functions for GFF manipulation.

def get_attr(attributes_string, key)
  attributes = attributes_string.split(';')
  attribute = attributes.grep(/#{key}=/)[0]
  attribute && attribute.split('=')[1]
end

def del_attr(attributes_string, key)
  attributes = attributes_string.split(';')
  attributes.delete_if do |attr|
    attr.match(/#{key}=/)
  end
  attributes.join(';')
end


